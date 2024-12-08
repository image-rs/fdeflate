use simd_adler32::Adler32;
use std::{
    collections::BinaryHeap,
    io::{self, Seek, SeekFrom, Write},
};

use crate::tables::{
    BITMASKS, CLCL_ORDER, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, LENGTH_TO_LEN_EXTRA,
    LENGTH_TO_SYMBOL,
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
    (0x330698ecu64.wrapping_mul(((v & 0xff_ffff) ^ 0x2722_0a95) as u64) >> 16) as u32
}
fn compute_hash4(v: u64) -> u32 {
    let mut hasher = fnv::FnvHasher::default();
    std::hash::Hasher::write_u32(&mut hasher, v as u32);
    std::hash::Hasher::finish(&hasher) as u32
}
fn compute_hash(v: u64) -> u32 {
    let mut hasher = fnv::FnvHasher::default();
    std::hash::Hasher::write_u64(&mut hasher, v);
    std::hash::Hasher::finish(&hasher) as u32
}

fn match_length(data: &[u8], anchor: usize, mut ip: usize, mut prev_index: usize) -> (u16, usize) {
    assert!(prev_index < ip, "{prev_index}, {ip}");

    if data[ip] != data[prev_index] {
        return (0, ip);
    }

    let mut length = 1;
    while length < 258 && ip > anchor && prev_index > 0 && data[ip - 1] == data[prev_index - 1] {
        length += 1;
        ip -= 1;
        prev_index -= 1;
    }
    while length < 258 && ip + length < data.len() && data[ip + length] == data[prev_index + length]
    {
        length += 1;
    }
    (length as u16, ip)
}

const CACHE3_SIZE: usize = 1 << 15;
const CACHE4_SIZE: usize = 1 << 15;
const CACHE5_SIZE: usize = 1 << 18;
const WINDOW_SIZE: usize = 32768;

struct CacheTable {
    hash3_table: Option<Box<[u32; CACHE3_SIZE]>>,
    hash4_table: Option<Box<[u32; CACHE4_SIZE]>>,
    hash_table: Box<[u32; CACHE5_SIZE]>,
    links: Box<[u32; WINDOW_SIZE]>,

    search_depth: u16,
    nice_length: u16,
    hash_mask: u64,
}
impl CacheTable {
    fn new(search_depth: u16, nice_length: u16, min_match: u8) -> Self {
        assert!((3..=6).contains(&min_match));

        Self {
            hash3_table: (min_match == 3)
                .then(|| vec![0; CACHE3_SIZE].into_boxed_slice().try_into().unwrap()),
            hash4_table: None, //Some(vec![0; CACHE4_SIZE].into_boxed_slice().try_into().unwrap()),
            hash_table: vec![0; CACHE5_SIZE].into_boxed_slice().try_into().unwrap(),
            links: vec![0; WINDOW_SIZE].into_boxed_slice().try_into().unwrap(),
            search_depth,
            nice_length,
            hash_mask: (1 << (min_match.max(4) * 8)) - 1,
        }
    }

    fn get_and_insert(
        &mut self,
        data: &[u8],
        anchor: usize,
        ip: usize,
        value: u64,
        min_match: u16,
    ) -> (u16, u16, usize) {
        let min_offset = ip.saturating_sub(32768).max(1);

        let mut best_offset = 0;
        let mut best_length = min_match - 1;
        let mut best_ip = 0;

        let mut n = self.search_depth;
        if min_match > 3 {
            n /= 2;
        }

        let hash = compute_hash(value & self.hash_mask);
        let hash_index = (hash as usize) % CACHE5_SIZE;
        let mut offset = self.hash_table[hash_index] as usize;

        // Insert current value
        self.hash_table[hash_index] = ip as u32;
        self.links[ip % WINDOW_SIZE] = offset as u32;

        // Visit previous matches
        loop {
            if offset < min_offset {
                break;
            }

            let (length, start) = match_length(data, anchor, ip, offset);
            if length > best_length {
                best_length = length;
                best_offset = offset as u32;
                best_ip = start;
            }
            if length >= self.nice_length || ip + length as usize == data.len() {
                break;
            }

            n -= 1;
            if n == 0 {
                break;
            }

            offset = self.links[offset % WINDOW_SIZE] as usize;
        }

        if let Some(hash4_table) = &mut self.hash4_table {
            let hash4 = compute_hash4(value);
            if best_length < 4 {
                let hash4_offset = hash4_table[(hash4 as usize) % CACHE4_SIZE] as usize;
                if hash4_offset >= ip.saturating_sub(8192).max(1) {
                    let (length, start) = match_length(data, anchor, ip, hash4_offset);
                    if length >= 4 {
                        best_length = length;
                        best_offset = hash4_offset as u32;
                        best_ip = start;
                    }
                }
            }
            hash4_table[(hash4 as usize) % CACHE4_SIZE] = ip as u32;

            if let Some(hash3_table) = &mut self.hash3_table {
                let hash3 = compute_hash3(value as u32);
                if best_length < min_match && min_match <= 3 {
                    let hash3_offset = hash3_table[(hash3 as usize) % CACHE3_SIZE] as usize;
                    if hash3_offset >= ip.saturating_sub(64).max(1) {
                        let (length, start) = match_length(data, anchor, ip, hash3_offset);
                        if length >= 3 {
                            best_length = length;
                            best_offset = hash3_offset as u32;
                            best_ip = start;
                        }
                    }
                }
                hash3_table[(hash3 as usize) % CACHE3_SIZE] = ip as u32;
            }
        }

        if best_length >= min_match {
            return (
                best_length as u16,
                (ip - best_offset as usize) as u16,
                best_ip,
            );
        }

        (0, 0, ip)
    }

    fn insert(&mut self, value: u64, offset: usize) {
        if let Some(hash4_table) = &mut self.hash4_table {
            let hash4 = compute_hash4(value);
            hash4_table[(hash4 as usize) % CACHE4_SIZE] = offset as u32;

            if let Some(hash3_table) = &mut self.hash3_table {
                let hash3 = compute_hash3(value as u32);
                hash3_table[(hash3 as usize) % CACHE3_SIZE] = offset as u32;
            }
        }

        let hash = compute_hash(value & self.hash_mask);
        let prev_offset = self.hash_table[(hash as usize) % CACHE5_SIZE];
        self.hash_table[(hash as usize) % CACHE5_SIZE] = offset as u32;
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

    min_match: u8,
    skip_ahead_shift: u8,
    search_depth: u16,
    nice_length: u16,
    max_lazy: u16,
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

            min_match: 4,
            search_depth: 1500,
            nice_length: 258,
            skip_ahead_shift: 29,
            max_lazy: 64,
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

        let mut matches = CacheTable::new(self.search_depth, self.nice_length, self.min_match);
        let mut ip = 0;

        while ip < data.len() {
            let mut symbols = Vec::new();

            let mut last_match = ip;
            'outer: while symbols.len() < 16384 && ip + 8 <= data.len() {
                let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());

                if current == 0 {
                    for j in last_match..ip {
                        symbols.push(Symbol::Literal(data[j]));
                    }

                    if ip == 0 || data[ip - 1] != 0 {
                        symbols.push(Symbol::Literal(0));
                        ip += 1;
                    }

                    let mut run_length = 0;
                    while ip < data.len() && data[ip] == 0 && run_length < 258 {
                        run_length += 1;
                        ip += 1;
                    }

                    symbols.push(Symbol::Backref {
                        length: run_length as u16,
                        distance: 1,
                        dist_sym: 0,
                    });

                    last_match = ip;
                    continue 'outer;
                }

                let (mut length, mut distance, mut match_start) =
                    matches.get_and_insert(&data, last_match, ip, current, 3);
                if length >= 3 {
                    if match_start == ip
                        && length < self.max_lazy
                        && ip + length as usize + 9 <= data.len()
                    {
                        ip += 1;
                        let (next_length, next_distance, next_match_start) =
                            matches.get_and_insert(&data, last_match, ip, current >> 8, length + 1);
                        if next_length > 0 {
                            if next_match_start < ip && distance > 1 {
                                distance = next_distance;
                                length = next_length;
                                match_start = next_match_start;
                            } else if next_match_start >= ip
                            /* && (next_length > length || (next_length == length && next_index * 4 < index)) */
                            {
                                length = next_length;
                                distance = next_distance;
                                match_start = next_match_start;
                            }
                        }
                    }
                    assert!(last_match <= match_start);

                    for j in last_match..match_start {
                        symbols.push(Symbol::Literal(data[j]));
                    }

                    symbols.push(Symbol::Backref {
                        length: length as u16,
                        distance,
                        dist_sym: distance_to_dist_sym(distance),
                    });

                    let match_end = match_start + length as usize;
                    let insert_end = (match_end - 2).min(data.len() - 8);
                    let insert_start = (ip + 1).max(insert_end.saturating_sub(16));
                    for j in (insert_start..insert_end).step_by(3) {
                        let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                        matches.insert(v, j);
                        matches.insert(v >> 8, j + 1);
                        matches.insert(v >> 16, j + 2);
                    }

                    // if insert_end >= insert_start + 3 {
                    //     let v = u64::from_le_bytes(data[insert_end - 3..][..8].try_into().unwrap());
                    //     matches.insert(v, insert_end - 3);
                    //     matches.insert(v >> 8, insert_end - 2);
                    //     matches.insert(v >> 16, insert_end - 1);
                    // }

                    ip = match_end;
                    last_match = ip;
                    continue 'outer;
                }

                // matches.insert(current >> 8, ip + 1);
                // matches.insert(current >> 16, ip + 2);

                // If we haven't found a match in a while, start skipping ahead by emitting multiple
                // literals at once.
                ip += 1 + ((ip - last_match) >> self.skip_ahead_shift);
            }
            if data.len() < ip + 8 {
                for j in last_match..data.len() {
                    symbols.push(Symbol::Literal(data[j]));
                }
                ip = data.len();
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

            let mut lengths = [0u8; 286];
            let mut codes = [0u16; 286];
            build_huffman_tree(&frequencies, &mut lengths, &mut codes, 15);

            let mut dist_lengths = [0u8; 30];
            let mut dist_codes = [0u16; 30];
            build_huffman_tree(&dist_frequencies, &mut dist_lengths, &mut dist_codes, 15);

            if ip == data.len() {
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
