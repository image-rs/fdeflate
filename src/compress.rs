use simd_adler32::Adler32;
use std::{
    collections::BinaryHeap,
    io::{self, Seek, SeekFrom, Write},
};

use crate::tables::{
    BITMASKS, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, HUFFMAN_CODES, HUFFMAN_LENGTHS,
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
    let mut dist_sym = 29;
    while dist_sym > 0 && distance < DIST_SYM_TO_DIST_BASE[dist_sym as usize - 1] {
        dist_sym -= 1;
    }
    dist_sym
}

fn hash(v: u64) -> u32 {
    (0x27220a95u64.wrapping_mul((v & 0xffff_ffff) ^ 0x56330698ec) >> 32) as u32
}
fn hash2(v: u64) -> u32 {
    (0x330698ecu64.wrapping_mul(v ^ 0x27220a95) >> 32) as u32
}

const WAYS: usize = 1;
const CACHE_SIZE: usize = 1 << 19;

#[derive(Debug, Copy, Clone)]
struct Entry {
    tags: [u32; WAYS],
    offsets: [u32; WAYS],
}

struct CacheTable {
    entries: Box<[Entry; CACHE_SIZE]>,
}
impl CacheTable {
    fn new() -> Self {
        let entries: Box<[Entry]> = vec![
            Entry {
                tags: [0; WAYS],
                offsets: [0; WAYS],
            };
            CACHE_SIZE
        ]
        .into_boxed_slice();

        Self {
            entries: entries.try_into().unwrap(),
        }
    }

    fn get(&self, data: &[u8], index: usize, hash: u32) -> (u32, u16) {
        let mut best_offset = 0;
        let mut best_length = 0;

        let entry = &self.entries[(hash as usize) % CACHE_SIZE];
        for i in 0..WAYS {
            if entry.tags[i] != hash {
                continue;
            }

            let offset = entry.offsets[i] as usize;
            if index - offset < 32768 {
                let mut length = 0;
                while length < 258
                    && index + length < data.len()
                    && data[index + length] == data[offset + length]
                {
                    length += 1;
                }

                if length > 3 && length > best_length {
                    best_offset = offset as u32;
                    best_length = length;
                }
                if length == 258 {
                    break;
                }
            }
        }

        (best_offset, best_length as u16)
    }

    fn contains(&self, hash: u32) -> bool {
        let entry = &self.entries[(hash as usize) % CACHE_SIZE];
        for i in 0..WAYS {
            if entry.tags[i] == hash {
                return true;
            }
        }
        false
    }

    fn insert(&mut self, hash: u32, offset: u32) {
        let entry = &mut self.entries[(hash as usize) % CACHE_SIZE];

        let mut oldest = 0;
        for i in 1..WAYS {
            if entry.offsets[i] < entry.offsets[oldest] {
                oldest = i;
            }
        }

        entry.tags[oldest] = hash;
        entry.offsets[oldest] = offset;
    }
}

/// Compressor that produces fdeflate compressed streams.
pub struct Compressor<W: Write> {
    checksum: Adler32,
    buffer: u64,
    nbits: u8,
    writer: W,
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
        self.writer.write_all(&HEADER[..53]).unwrap();
        self.write_bits(HEADER[53] as u64, 5)?;

        Ok(())
    }

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);

        enum Symbol {
            Literal(u8),
            Rle {
                length: u16,
            },
            Backref {
                length: u16,
                distance: u16,
                dist_sym: u8,
            },
        }

        let mut short_matches = CacheTable::new();
        let mut long_matches = CacheTable::new();
        let mut i = 0;

        let mut lengths = HUFFMAN_LENGTHS;
        let mut dist_lengths = [6u8; 30];
        dist_lengths[0] = 1;

        while i < data.len() {
            let mut symbols = Vec::new();

            let block_end = data.len().min(i + 128 * 1024);

            for len in &mut lengths {
                if *len == 0 { *len = 15; }
            }
            for len in &mut dist_lengths {
                if *len == 0 { *len = 6; }
            }

            let mut last_match = i;
            while i < block_end && i + 8 < data.len() {
                let current = u64::from_le_bytes(data[i..][..8].try_into().unwrap());

                if current & 0xffff_ffff == 0 {
                    let mut run_length = 4;
                    while run_length <= 258
                        && i + run_length < data.len()
                        && data[i + run_length] == 0
                    {
                        run_length += 1;
                    }
                    symbols.push(Symbol::Literal(0));
                    symbols.push(Symbol::Rle {
                        length: (run_length - 1) as u16,
                    });
                    i += run_length;
                    continue;
                }

                // Long hash
                let long_hash = hash2(current);
                let (prev_i, length) = long_matches.get(data, i, long_hash as u32);
                if length >= 8 {
                    for j in (i + length as usize - 8)..(i + length as usize).min(data.len() - 8) {
                        let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                        short_matches.insert(hash(v), j as u32);
                        long_matches.insert(hash2(v), j as u32);
                    }

                    let dist = (i - prev_i as usize) as u16;
                    let dist_sym = distance_to_dist_sym(dist);

                    symbols.push(Symbol::Backref {
                        length: length as u16,
                        distance: dist,
                        dist_sym,
                    });
                    i += length as usize;
                    last_match = i;
                    continue;
                }

                // Short hash
                let short_hash = hash(current);
                let (prev_i, length) = short_matches.get(data, i, short_hash as u32);
                short_matches.insert(short_hash as u32, i as u32);
                long_matches.insert(long_hash as u32, i as u32);
                if length >= 3 {
                    let next = u64::from_le_bytes(data[i + 1..][..8].try_into().unwrap());
                    if next != 0 && !long_matches.contains(hash2(next)) {
                        let dist = (i - prev_i as usize) as u16;
                        let dist_sym = distance_to_dist_sym(dist);

                        let sym = LENGTH_TO_SYMBOL[length as usize - 3] as usize;
                        let len_bits = lengths[sym];
                        let len_extra = LENGTH_TO_LEN_EXTRA[length as usize - 3];
                        let dist_bits = dist_lengths[dist_sym as usize];
                        let dist_extra = DIST_SYM_TO_DIST_EXTRA[dist_sym as usize];
                        let backref_cost = (len_bits + len_extra + dist_bits + dist_extra) as u32;

                        let mut literal_cost = 0;
                        for j in i..i + length as usize {
                            literal_cost += lengths[data[j] as usize] as u32;
                            if literal_cost >= backref_cost {
                                break;
                            }
                        }
                        if literal_cost > backref_cost {
                            symbols.push(Symbol::Backref {
                                length: length as u16,
                                distance: dist,
                                dist_sym,
                            });

                            for j in (i + 1).min(i + length as usize).saturating_sub(8)
                                ..(i + length as usize).min(data.len() - 8)
                            {
                                let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                                short_matches.insert(hash(v), j as u32);
                                long_matches.insert(hash2(v), j as u32);
                            }

                            i += length as usize;
                            last_match = i;
                            continue;
                        }
                    }
                }

                for _ in 0..((i - last_match) >> 9).min((data.len() - i).saturating_sub(8)) {
                    symbols.push(Symbol::Literal(data[i]));
                    i += 1;
                }

                symbols.push(Symbol::Literal(data[i]));
                i += 1;
            }
            for i in i..block_end {
                symbols.push(Symbol::Literal(data[i]));
            }
            i = i.max(block_end);

            let mut frequencies = [0u32; 286];
            let mut dist_frequencies = [0u32; 30];
            frequencies[256] = 1;
            for symbol in &symbols {
                match symbol {
                    Symbol::Literal(lit) => frequencies[*lit as usize] += 1,
                    Symbol::Rle { length } => {
                        let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                        frequencies[sym] += 1;
                        dist_frequencies[0] += 1;
                    }
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

            let rle_three = lengths[257] + dist_lengths[0] > 3 * lengths[0];
            let rle_four = lengths[257] + dist_lengths[0] > 4 * lengths[0];

            for symbol in &symbols {
                match symbol {
                    Symbol::Literal(lit) => {
                        let sym = *lit as usize;
                        self.write_bits(codes[sym] as u64, lengths[sym] as u8)?;
                    }
                    Symbol::Rle { length: 3 } if rle_three => {
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                    }
                    Symbol::Rle { length: 4 } if rle_four => {
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                        self.write_bits(codes[0] as u64, lengths[0] as u8)?;
                    }
                    Symbol::Rle { length } => {
                        let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                        self.write_bits(codes[sym] as u64, lengths[sym] as u8)?;
                        let len_extra = LENGTH_TO_LEN_EXTRA[*length as usize - 3];
                        let extra = (((*length as u32) - 3) & BITMASKS[len_extra as usize]) as u64;
                        self.write_bits(extra, len_extra + 1)?;
                        self.write_bits(dist_codes[0] as u64, dist_lengths[0])?;
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
                        let extra = ((*distance as u32) & BITMASKS[dist_extra as usize]) as u64;
                        self.write_bits(extra, dist_extra)?;
                    }
                }
            }
            self.write_bits(codes[256] as u64, lengths[256])?;
        }

        Ok(())
    }

    /// Write the remainder of the stream and return the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        // Write end of block
        self.write_bits(HUFFMAN_CODES[256] as u64, HUFFMAN_LENGTHS[256])?;
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
    use super::*;
    use rand::Rng;

    fn roundtrip(data: &[u8]) {
        let compressed = compress_to_vec(data);
        let decompressed = miniz_oxide::inflate::decompress_to_vec_zlib(&compressed).unwrap();
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
