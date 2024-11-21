use simd_adler32::Adler32;
use std::{
    collections::BinaryHeap,
    io::{self, Seek, SeekFrom, Write},
};

use crate::tables::{
    BITMASKS, CLCL_ORDER, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, HUFFMAN_LENGTHS,
    LENGTH_TO_LEN_EXTRA, LENGTH_TO_SYMBOL,
};

fn build_huffman_tree(
    frequencies: &[u32],
    lengths: &mut [u8],
    codes: &mut [u16],
    length_limit: u8,
) -> bool {
    assert_eq!(frequencies.len(), lengths.len());
    assert_eq!(frequencies.len(), codes.len());

    if frequencies.iter().filter(|&&f| f > 0).count() <= 1 {
        lengths.fill(0);
        codes.fill(0);
        if let Some(i) = frequencies.iter().position(|&f| f > 0) {
            lengths[i] = 1;
        }
        return false;
    }

    #[derive(Eq, PartialEq, Copy, Clone, Debug)]
    struct Item(u32, u16);
    impl Ord for Item {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            other.0.cmp(&self.0)
        }
    }
    impl PartialOrd for Item {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    // Build a huffman tree
    let mut internal_nodes = Vec::new();
    let mut nodes = BinaryHeap::from_iter(
        frequencies
            .iter()
            .enumerate()
            .filter(|(_, &frequency)| frequency > 0)
            .map(|(i, &frequency)| Item(frequency, i as u16)),
    );
    while nodes.len() > 1 {
        let Item(frequency1, index1) = nodes.pop().unwrap();
        let mut root = nodes.peek_mut().unwrap();
        internal_nodes.push((index1, root.1));
        *root = Item(
            frequency1 + root.0,
            internal_nodes.len() as u16 + frequencies.len() as u16 - 1,
        );
    }

    // Walk the tree to assign code lengths
    lengths.fill(0);
    let mut stack = Vec::new();
    stack.push((nodes.pop().unwrap().1, 0));
    while let Some((node, depth)) = stack.pop() {
        let node = node as usize;
        if node < frequencies.len() {
            lengths[node] = depth as u8;
        } else {
            let (left, right) = internal_nodes[node - frequencies.len()];
            stack.push((left, depth + 1));
            stack.push((right, depth + 1));
        }
    }

    // Limit the codes to length length_limit
    let mut max_length = 0;
    for &length in lengths.iter() {
        max_length = max_length.max(length);
    }
    if max_length > length_limit {
        let mut counts = [0u32; 16];
        for &length in lengths.iter() {
            counts[length.min(length_limit) as usize] += 1;
        }

        let mut total = 0;
        for (i, count) in counts
            .iter()
            .enumerate()
            .skip(1)
            .take(length_limit as usize)
        {
            total += count << (length_limit as usize - i);
        }

        while total > 1u32 << length_limit {
            let mut i = length_limit as usize - 1;
            while counts[i] == 0 {
                i -= 1;
            }
            counts[i] -= 1;
            counts[length_limit as usize] -= 1;
            counts[i + 1] += 2;
            total -= 1;
        }

        // assign new lengths
        let mut len = length_limit;
        let mut indexes = frequencies.iter().copied().enumerate().collect::<Vec<_>>();
        indexes.sort_unstable_by_key(|&(_, frequency)| frequency);
        for &(i, frequency) in indexes.iter() {
            if frequency > 0 {
                while counts[len as usize] == 0 {
                    len -= 1;
                }
                lengths[i] = len;
                counts[len as usize] -= 1;
            }
        }
    }

    // Assign codes
    codes.fill(0);
    let mut code = 0u32;
    for len in 1..=length_limit {
        for (i, &length) in lengths.iter().enumerate() {
            if length == len {
                codes[i] = (code as u16).reverse_bits() >> (16 - len);
                code += 1;
            }
        }
        code <<= 1;
    }
    assert_eq!(code, 2 << length_limit);

    true
}

fn distance_to_dist_sym(distance: u16) -> u8 {
    const LOOKUP: [u8; 16] = [0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7];
    if distance <= 16 {
        return LOOKUP[distance as usize - 1];
    }

    let mut dist_sym = 29;
    while dist_sym > 0 && distance < DIST_SYM_TO_DIST_BASE[dist_sym as usize] {
        dist_sym -= 1;
    }
    dist_sym
}

fn compute_hash3(v: u32) -> u32 {
    (0x330698ecu64.wrapping_mul(((v & 0xff_ffff) ^ 0x2722_0a95) as u64) >> 32) as u32
}
fn compute_hash4(v: u32) -> u32 {
    (0x27220a95u64.wrapping_mul((v ^ 0x3306_98ec) as u64) >> 32) as u32
}

fn match_length(data: &[u8], index: usize, prev_index: usize) -> u16 {
    assert!(prev_index < index);

    let mut length = 0;
    while length < 258
        && index + length < data.len()
        && data[index + length] == data[prev_index + length]
    {
        length += 1;
    }
    length as u16
}

const CACHE3_SIZE: usize = 1 << 14;
const CACHE4_SIZE: usize = 1 << 16;
const WINDOW_SIZE: usize = 32768;

const SEARCH_DEPTH: u16 = 100;
const NICE_LENGTH: u16 = 60;

struct CacheTable {
    hash3_table: Box<[u32; CACHE3_SIZE]>,
    hash4_table: Box<[u32; CACHE4_SIZE]>,
    links: Box<[u32; WINDOW_SIZE]>,
}
impl CacheTable {
    fn new() -> Self {
        Self {
            hash3_table: vec![0; CACHE3_SIZE].into_boxed_slice().try_into().unwrap(),
            hash4_table: vec![0; CACHE4_SIZE].into_boxed_slice().try_into().unwrap(),
            links: vec![0; WINDOW_SIZE].into_boxed_slice().try_into().unwrap(),
        }
    }

    fn get(&self, data: &[u8], index: usize, hash3: u32, hash4: u32, min_match: u16) -> (u32, u16) {
        // if index + min_match as usize >= data.len() {
        //     return (0, 0);
        // }

        let min_offset = index.saturating_sub(32768).max(1);

        let mut best_offset = 0;
        let mut best_length = min_match - 1;

        // if min_match == 3 {
        //     let hash3_offset = self.hash3_table[(hash3 as usize) % CACHE3_SIZE] as usize;
        //     if hash3_offset >= index.saturating_sub(8192).max(1) {
        //         best_length = match_length(data, index, hash3_offset);
        //         best_offset = hash3_offset as u32;
        //     }
        // }

        let mut n = SEARCH_DEPTH;
        let mut offset = self.hash4_table[(hash4 as usize) % CACHE4_SIZE] as usize;
        loop {
            if offset < min_offset {
                break;
            }

            // if data[index + best_length as usize] == data[offset + best_length as usize] {
            let length = match_length(data, index, offset);
            if length > best_length {
                best_length = length;
                best_offset = offset as u32;
            }
            if length >= NICE_LENGTH || index + length as usize == data.len() {
                break;
            }
            // }

            n -= 1;
            if n == 0 {
                break;
            }

            offset = self.links[offset % WINDOW_SIZE] as usize;
        }

        if best_length >= min_match {
            return (best_offset as u32, best_length as u16);
        }

        (0, 0)
    }

    fn insert(&mut self, hash3: u32, hash4: u32, offset: usize) {
        // self.hash3_table[(hash3 as usize) % CACHE3_SIZE] = offset as u32;

        let prev_offset = self.hash4_table[(hash4 as usize) % CACHE4_SIZE];
        self.hash4_table[(hash4 as usize) % CACHE4_SIZE] = offset as u32;
        self.links[offset as usize % WINDOW_SIZE] = prev_offset;
    }
}

/// Compressor that produces fdeflate compressed streams.
pub struct Compressor<W: Write> {
    checksum: Adler32,
    buffer: u64,
    nbits: u8,
    writer: W,
    pending: Vec<u8>,
}
impl<W: Write> Compressor<W> {
    fn write_bits(&mut self, bits: u64, nbits: u8) -> io::Result<()> {
        debug_assert!(nbits <= 64);

        self.buffer |= bits << self.nbits;
        self.nbits += nbits;

        if self.nbits >= 64 {
            self.writer.write_all(&self.buffer.to_le_bytes())?;
            self.nbits -= 64;
            self.buffer = bits.checked_shr((nbits - self.nbits) as u32).unwrap_or(0);
        }
        debug_assert!(self.nbits < 64);
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.nbits % 8 != 0 {
            self.write_bits(0, 8 - self.nbits % 8)?;
        }
        if self.nbits > 0 {
            self.writer
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])
                .unwrap();
            self.buffer = 0;
            self.nbits = 0;
        }
        Ok(())
    }

    /// Create a new Compressor.
    pub fn new(writer: W) -> io::Result<Self> {
        let mut compressor = Self {
            checksum: Adler32::new(),
            buffer: 0,
            nbits: 0,
            writer,
            pending: Vec::new(),
        };
        compressor.write_headers()?;
        Ok(compressor)
    }

    fn write_headers(&mut self) -> io::Result<()> {
        const HEADER: [u8; 54] = [
            120, 1, 237, 192, 3, 160, 36, 89, 150, 198, 241, 255, 119, 238, 141, 200, 204, 167,
            114, 75, 99, 174, 109, 219, 182, 109, 219, 182, 109, 219, 182, 109, 105, 140, 158, 150,
            74, 175, 158, 50, 51, 34, 238, 249, 118, 183, 106, 122, 166, 135, 59, 107, 213, 15,
        ];
        self.writer.write_all(&HEADER[..2]).unwrap();
        //self.write_bits(HEADER[53] as u64, 5)?;

        Ok(())
    }

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);
        self.pending.extend_from_slice(data);
        Ok(())
    }

    /// Write the remainder of the stream and return the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        let data = std::mem::take(&mut self.pending);

        enum Symbol {
            Literal(u8),
            Backref {
                length: u16,
                distance: u16,
                dist_sym: u8,
            },
        }

        let mut matches = CacheTable::new();
        let mut i = 0;

        let mut lengths = HUFFMAN_LENGTHS;
        let mut dist_lengths = [5u8; 30];
        //dist_lengths[0] = 1;

        while i < data.len() {
            let mut symbols = Vec::new();

            // let block_end = data.len().min(i + 64 * 1024);

            for len in &mut lengths {
                if *len == 0 {
                    *len = 15;
                }
            }
            for len in &mut dist_lengths {
                if *len == 0 {
                    *len = 5;
                }
            }

            let mut last_match = i;
            'outer: while symbols.len() < 16384 && i + 4 <= data.len() {
                let current = u32::from_le_bytes(data[i..][..4].try_into().unwrap());

                let hash3 = compute_hash3(current);
                let hash4 = compute_hash4(current);
                let (mut index, mut length) = matches.get(&data, i, hash3, hash4, 3);
                matches.insert(hash3, hash4, i);

                if length >= 3 {
                    'peak_ahead: loop {
                        if length < NICE_LENGTH && i + length as usize + 4 <= data.len() {
                            let next = u32::from_le_bytes(data[i + 1..][..4].try_into().unwrap());
                            let next_hash4 = compute_hash4(next);
                            let (next_index, next_length) =
                                matches.get(&data, i + 1, 0, next_hash4, length + 1);

                            if next_length > length || (next_length == length && next_index * 4 < index) {
                                // let next2 =
                                //     u32::from_le_bytes(data[i + 2..][..4].try_into().unwrap());
                                // let next2_hash4 = compute_hash4(next);
                                // let (next2_index, next2_length) =
                                //     matches.get(&data, i + 2, 0, next2_hash4, length + 1);

                                // if next2_length > next_length {
                                //     let next3 =
                                //         u32::from_le_bytes(data[i + 3..][..4].try_into().unwrap());
                                //     let next3_hash4 = compute_hash4(next);
                                //     let (next3_index, next3_length) =
                                //         matches.get(&data, i + 3, 0, next3_hash4, length + 1);

                                //     if next3_length > next2_length {
                                //         let next4 = u32::from_le_bytes(
                                //             data[i + 4..][..4].try_into().unwrap(),
                                //         );
                                //         let next4_hash4 = compute_hash4(next);
                                //         let (next4_index, next4_length) =
                                //             matches.get(&data, i + 4, 0, next4_hash4, length + 1);

                                //         if next4_length > next3_length {
                                //             let dist = (i - index as usize) as u16;
                                //             symbols.push(Symbol::Backref {
                                //                 length: 4,
                                //                 distance: dist,
                                //                 dist_sym: distance_to_dist_sym(dist),
                                //             });

                                //             matches.insert(compute_hash3(next), next_hash4, i);
                                //             matches.insert(
                                //                 compute_hash3(next2),
                                //                 next2_hash4,
                                //                 i + 1,
                                //             );
                                //             matches.insert(
                                //                 compute_hash3(next3),
                                //                 next3_hash4,
                                //                 i + 2,
                                //             );
                                //             matches.insert(
                                //                 compute_hash3(next4),
                                //                 next4_hash4,
                                //                 i + 3,
                                //             );

                                //             i += 3;
                                //             last_match = i;
                                //             index = next4_index;
                                //             length = next4_length;
                                //             continue 'peak_ahead;
                                //         } else {
                                //             symbols.push(Symbol::Literal(data[i]));
                                //             symbols.push(Symbol::Literal(data[i + 1]));
                                //             symbols.push(Symbol::Literal(data[i + 2]));
                                //             matches.insert(compute_hash3(next), next_hash4, i);
                                //             matches.insert(compute_hash3(next2), next2_hash4, i + 1);
                                //             matches.insert(compute_hash3(next3), next3_hash4, i + 2);

                                //             i += 3;
                                //             index = next3_index;
                                //             length = next3_length;
                                //         }
                                //     } else {
                                //         symbols.push(Symbol::Literal(data[i]));
                                //         symbols.push(Symbol::Literal(data[i + 1]));
                                //         matches.insert(compute_hash3(next), next_hash4, i);
                                //         matches.insert(compute_hash3(next2), next2_hash4, i + 1);

                                //         i += 2;
                                //         index = next2_index;
                                //         length = next2_length;
                                //     }
                                // } else {
                                symbols.push(Symbol::Literal(data[i]));
                                matches.insert(compute_hash3(next), next_hash4, i);

                                i += 1;
                                index = next_index;
                                length = next_length;
                                // }
                            }
                        }

                        let dist = (i - index as usize) as u16;
                        symbols.push(Symbol::Backref {
                            length: length as u16,
                            distance: dist,
                            dist_sym: distance_to_dist_sym(dist),
                        });

                        for j in (i + 1)..(i + length as usize).min(data.len() - 4) {
                            let v = u32::from_le_bytes(data[j..][..4].try_into().unwrap());
                            matches.insert(compute_hash3(v), compute_hash4(v), j);
                        }

                        i += length as usize;
                        last_match = i;
                        continue 'outer;
                    }
                }

                // // If we haven't found a match in a while, start skipping ahead by emitting multiple
                // // literals at once.
                // for _ in 0..((i - last_match) >> 9).min((data.len() - i).saturating_sub(1)) {
                //     symbols.push(Symbol::Literal(data[i]));
                //     i += 1;
                // }

                // Emit a literal
                symbols.push(Symbol::Literal(data[i]));
                i += 1;
            }
            if data.len() - i < 4 {
                for j in i..data.len() {
                    symbols.push(Symbol::Literal(data[j]));
                }
                i = data.len();
            }

            let mut frequencies = [0u32; 286];
            let mut dist_frequencies = [0u32; 30];
            frequencies[256] = 1;
            for symbol in &symbols {
                match symbol {
                    Symbol::Literal(lit) => frequencies[*lit as usize] += 1,
                    Symbol::Backref {
                        length, dist_sym, ..
                    } => {
                        let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                        frequencies[sym] += 1;
                        dist_frequencies[*dist_sym as usize] += 1;
                    }
                }
            }

            let mut codes = [0u16; 286];
            build_huffman_tree(&frequencies, &mut lengths, &mut codes, 15);

            let mut dist_codes = [0u16; 30];
            build_huffman_tree(&dist_frequencies, &mut dist_lengths, &mut dist_codes, 15);

            if i == data.len() {
                self.write_bits(101, 3)?; // final block
            } else {
                self.write_bits(100, 3)?; // non-final block
            }
            self.write_bits(29, 5)?; // hlit
            self.write_bits(29, 5)?; // hdist
            self.write_bits(15, 4)?; // hclen

            let mut code_length_frequencies = [0u32; 19];
            for &length in &lengths {
                code_length_frequencies[length as usize] += 1;
            }
            for &length in &dist_lengths {
                code_length_frequencies[length as usize] += 1;
            }
            let mut code_length_lengths = [0u8; 19];
            let mut code_length_codes = [0u16; 19];
            build_huffman_tree(
                &code_length_frequencies,
                &mut code_length_lengths,
                &mut code_length_codes,
                7,
            );

            for j in 0..19 {
                self.write_bits(code_length_lengths[CLCL_ORDER[j]] as u64, 3)?;
            }

            for &length in lengths.iter().chain(&dist_lengths) {
                self.write_bits(
                    code_length_codes[length as usize] as u64,
                    code_length_lengths[length as usize],
                )?;
            }

            for symbol in &symbols {
                match symbol {
                    Symbol::Literal(lit) => {
                        let sym = *lit as usize;
                        self.write_bits(codes[sym] as u64, lengths[sym] as u8)?;
                    }
                    Symbol::Backref {
                        length,
                        distance,
                        dist_sym,
                    } => {
                        let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                        self.write_bits(codes[sym] as u64, lengths[sym] as u8)?;
                        let len_extra = LENGTH_TO_LEN_EXTRA[*length as usize - 3];
                        let extra = (((*length as u32) - 3) & BITMASKS[len_extra as usize]) as u64;
                        self.write_bits(extra, len_extra)?;

                        self.write_bits(
                            dist_codes[*dist_sym as usize] as u64,
                            dist_lengths[*dist_sym as usize],
                        )?;
                        let dist_extra = DIST_SYM_TO_DIST_EXTRA[*dist_sym as usize];
                        let extra = *distance - DIST_SYM_TO_DIST_BASE[*dist_sym as usize];

                        self.write_bits(extra as u64, dist_extra)?;
                    }
                }
            }
            self.write_bits(codes[256] as u64, lengths[256])?;
        }

        // Write end of block
        self.flush()?;

        // Write Adler32 checksum
        let checksum: u32 = self.checksum.finish();
        self.writer
            .write_all(checksum.to_be_bytes().as_ref())
            .unwrap();
        Ok(self.writer)
    }
}

/// Compressor that only writes the stored blocks.
///
/// This is useful for writing files that are not compressed, but still need to be wrapped in a
/// zlib stream.
pub struct StoredOnlyCompressor<W> {
    writer: W,
    checksum: Adler32,
    block_bytes: u16,
}
impl<W: Write + Seek> StoredOnlyCompressor<W> {
    /// Creates a new `StoredOnlyCompressor` that writes to the given writer.
    pub fn new(mut writer: W) -> io::Result<Self> {
        writer.write_all(&[0x78, 0x01])?; // zlib header
        writer.write_all(&[0; 5])?; // placeholder stored block header

        Ok(Self {
            writer,
            checksum: Adler32::new(),
            block_bytes: 0,
        })
    }

    fn set_block_header(&mut self, size: u16, last: bool) -> io::Result<()> {
        self.writer.seek(SeekFrom::Current(-(size as i64 + 5)))?;
        self.writer.write_all(&[
            last as u8,
            (size & 0xFF) as u8,
            ((size >> 8) & 0xFF) as u8,
            (!size & 0xFF) as u8,
            ((!size >> 8) & 0xFF) as u8,
        ])?;
        self.writer.seek(SeekFrom::Current(size as i64))?;

        Ok(())
    }

    /// Writes the given data to the underlying writer.
    pub fn write_data(&mut self, mut data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);
        while !data.is_empty() {
            if self.block_bytes == u16::MAX {
                self.set_block_header(u16::MAX, false)?;
                self.writer.write_all(&[0; 5])?; // placeholder stored block header
                self.block_bytes = 0;
            }

            let prefix_bytes = data.len().min((u16::MAX - self.block_bytes) as usize);
            self.writer.write_all(&data[..prefix_bytes])?;
            self.block_bytes += prefix_bytes as u16;
            data = &data[prefix_bytes..];
        }

        Ok(())
    }

    /// Finish writing the final block and return the underlying writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.set_block_header(self.block_bytes, true)?;

        // Write Adler32 checksum
        let checksum: u32 = self.checksum.finish();
        self.writer
            .write_all(checksum.to_be_bytes().as_ref())
            .unwrap();

        Ok(self.writer)
    }
}
impl<W> StoredOnlyCompressor<W> {
    /// Return the number of bytes that will be written to the output stream
    /// for the given input size. Because this compressor only writes stored blocks,
    /// the output size is always slightly *larger* than the input size.
    pub fn compressed_size(raw_size: usize) -> usize {
        (raw_size.saturating_sub(1) / u16::MAX as usize) * (u16::MAX as usize + 5)
            + (raw_size % u16::MAX as usize + 5)
            + 6
    }
}

/// Compresses the given data.
pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new(Vec::with_capacity(input.len() / 4)).unwrap();
    compressor.write_data(input).unwrap();
    compressor.finish().unwrap()
}

#[cfg(test)]
mod tests {
    use crate::decompress;

    use super::*;
    use rand::Rng;

    #[test]
    fn test_distance_to_dist_sym() {
        assert_eq!(distance_to_dist_sym(1), 0);
        assert_eq!(distance_to_dist_sym(2), 1);
        assert_eq!(distance_to_dist_sym(3), 2);
        assert_eq!(distance_to_dist_sym(4), 3);
        assert_eq!(distance_to_dist_sym(5), 4);
        assert_eq!(distance_to_dist_sym(7), 5);
        assert_eq!(distance_to_dist_sym(9), 6);
        assert_eq!(distance_to_dist_sym(13), 7);
        assert_eq!(distance_to_dist_sym(18), 8);
        assert_eq!(distance_to_dist_sym(257), 16);
    }

    fn roundtrip(data: &[u8]) {
        let compressed = compress_to_vec(data);
        //let decompressed = miniz_oxide::inflate::decompress_to_vec_zlib(&compressed).unwrap();
        let decompressed = crate::decompress_to_vec(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn it_works() {
        roundtrip(b"Hello world!");
    }

    #[test]
    fn constant() {
        roundtrip(&vec![0; 2048]);
        roundtrip(&vec![5; 2048]);
        roundtrip(&vec![128; 2048]);
        roundtrip(&vec![254; 2048]);
    }

    #[test]
    fn random() {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 2048];
        for _ in 0..10 {
            for byte in &mut data {
                *byte = rng.gen();
            }
            roundtrip(&data);
        }
    }
}
