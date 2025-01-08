use hc_matchfinder::HashChainMatchFinder;
use simd_adler32::Adler32;
use std::{
    collections::BinaryHeap,
    io::{self, Seek, SeekFrom, Write},
};

use crate::tables::{
    BITMASKS, CLCL_ORDER, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, LENGTH_TO_LEN_EXTRA,
    LENGTH_TO_SYMBOL,
};

mod hc_matchfinder;

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
fn compute_hash(v: u64) -> u32 {
    let mut hasher = fnv::FnvHasher::default();
    std::hash::Hasher::write_u64(&mut hasher, v);
    std::hash::Hasher::finish(&hasher) as u32
}

enum Symbol {
    LiteralRun {
        start: u32,
        end: u32,
    },
    Backref {
        length: u16,
        distance: u16,
        dist_sym: u8,
    },
}

enum CompressorInner {
    Stored,
    Fast {
        min_match: u8,
        skip_ahead_shift: u8,
        search_depth: u16,
        nice_length: u16,
    },
    Slow {
        search_depth: u16,
        nice_length: u16,
        max_lazy: u16,
    },
}
impl CompressorInner {
    fn compress_data<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        eof: bool,
    ) -> io::Result<()> {
        match self {
            Self::Stored => Self::compress_stored(writer, data, eof),
            Self::Fast {
                min_match,
                skip_ahead_shift,
                search_depth,
                nice_length,
            } => Self::compress_fast(
                writer,
                data,
                *min_match,
                *skip_ahead_shift,
                *search_depth,
                *nice_length,
            ),
            Self::Slow {
                search_depth,
                nice_length,
                max_lazy,
            } => {
                todo!()
            }
        }
    }

    fn compress_stored<W: Write>(
        writer: &mut BitWriter<W>,
        data: &[u8],
        eof: bool,
    ) -> io::Result<()> {
        if data.is_empty() {
            if eof {
                // TODO: write empty final block
            }
            return Ok(());
        }

        let chunks = data.chunks(65535);
        let last_chunk_index = chunks.len() - 1;
        for (i, chunk) in chunks.into_iter().enumerate() {
            if i == last_chunk_index {
                writer.write_bits(1, 3)?; // final block
            } else {
                writer.write_bits(0, 3)?; // non-final block
            }
            writer.flush()?;
            writer
                .writer
                .write_all(&(chunk.len() as u16).to_le_bytes())?;
            writer
                .writer
                .write_all(&(!(chunk.len() as u16)).to_le_bytes())?;
            writer.writer.write_all(chunk)?;
        }
        return Ok(());
    }

    fn write_symbols<W: Write>(
        writer: &mut BitWriter<W>,
        data: &[u8],
        symbols: &[Symbol],
        eof: bool,
    ) -> io::Result<()> {
        let mut frequencies = [0u32; 286];
        let mut dist_frequencies = [0u32; 30];
        frequencies[256] = 1;
        for symbol in symbols {
            match symbol {
                Symbol::LiteralRun { start, end } => {
                    for lit in &data[*start as usize..*end as usize] {
                        frequencies[*lit as usize] += 1;
                    }
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

        let mut lengths = [0u8; 286];
        let mut codes = [0u16; 286];
        build_huffman_tree(&frequencies, &mut lengths, &mut codes, 15);

        let mut dist_lengths = [0u8; 30];
        let mut dist_codes = [0u16; 30];
        build_huffman_tree(&dist_frequencies, &mut dist_lengths, &mut dist_codes, 15);

        if eof {
            writer.write_bits(101, 3)?; // final block
        } else {
            writer.write_bits(100, 3)?; // non-final block
        }
        writer.write_bits(29, 5)?; // hlit
        writer.write_bits(29, 5)?; // hdist
        writer.write_bits(15, 4)?; // hclen

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
            writer.write_bits(code_length_lengths[CLCL_ORDER[j]] as u64, 3)?;
        }

        for &length in lengths.iter().chain(&dist_lengths) {
            writer.write_bits(
                code_length_codes[length as usize] as u64,
                code_length_lengths[length as usize],
            )?;
        }

        for symbol in symbols {
            match symbol {
                Symbol::LiteralRun { start, end } => {
                    let mut groups = data[*start as usize..*end as usize].chunks_exact(4);
                    for group in &mut groups {
                        let code0 = codes[group[0] as usize] as u64;
                        let code1 = codes[group[1] as usize] as u64;
                        let code2 = codes[group[2] as usize] as u64;
                        let code3 = codes[group[3] as usize] as u64;

                        let len0 = lengths[group[0] as usize];
                        let len1 = lengths[group[1] as usize];
                        let len2 = lengths[group[2] as usize];
                        let len3 = lengths[group[3] as usize];

                        writer.write_bits(
                            code0
                                | (code1 << len0)
                                | (code2 << (len0 + len1))
                                | (code3 << (len0 + len1 + len2)),
                            len0 + len1 + len2 + len3,
                        )?;
                    }

                    for &lit in groups.remainder() {
                        writer
                            .write_bits(codes[lit as usize] as u64, lengths[lit as usize] as u8)?;
                    }
                }
                Symbol::Backref {
                    length,
                    distance,
                    dist_sym,
                } => {
                    let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                    writer.write_bits(codes[sym] as u64, lengths[sym] as u8)?;
                    let len_extra = LENGTH_TO_LEN_EXTRA[*length as usize - 3];
                    let extra = (((*length as u32) - 3) & BITMASKS[len_extra as usize]) as u64;
                    writer.write_bits(extra, len_extra)?;

                    writer.write_bits(
                        dist_codes[*dist_sym as usize] as u64,
                        dist_lengths[*dist_sym as usize],
                    )?;
                    let dist_extra = DIST_SYM_TO_DIST_EXTRA[*dist_sym as usize];
                    let extra = *distance - DIST_SYM_TO_DIST_BASE[*dist_sym as usize];

                    writer.write_bits(extra as u64, dist_extra)?;
                }
            }
        }
        writer.write_bits(codes[256] as u64, lengths[256])?;
        Ok(())
    }

    fn compress_fast<W: Write>(
        writer: &mut BitWriter<W>,
        data: &[u8],
        min_match: u8,
        skip_ahead_shift: u8,
        search_depth: u16,
        nice_length: u16,
    ) -> io::Result<()> {
        let mut matches = HashChainMatchFinder::new(search_depth, nice_length, min_match);
        let mut ip = 0;

        let mut length = 0;
        let mut distance = 0;
        let mut match_start = 0;

        while ip < data.len() {
            let mut symbols = Vec::new();

            let mut last_match = ip;
            'outer: while symbols.len() < 16384 && ip + 8 <= data.len() {
                let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());

                if length == 0 {
                    if current == 0 {
                        while ip > last_match && data[ip - 1] == 0 {
                            ip -= 1;
                        }

                        if ip == 0 || data[ip - 1] != 0 {
                            ip += 1;
                        }

                        symbols.push(Symbol::LiteralRun {
                            start: last_match as u32,
                            end: ip as u32,
                        });

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

                        length = 0;
                        continue;
                    }

                    (length, distance, match_start) =
                        matches.get_and_insert(&data, last_match, ip, current, 3);
                }

                if length >= 3 {
                    // if match_start + length as usize > ip + 1
                    //     && length < max_lazy
                    //     && ip + length as usize + 9 <= data.len()
                    // {
                    //     ip += 1;
                    //     let (next_length, next_distance, next_match_start) =
                    //         matches.get_and_insert(&data, last_match, ip, current >> 8, length + 1);
                    //     if next_length > 0 && match_start + 1 >= next_match_start {
                    //         distance = next_distance;
                    //         length = next_length;
                    //         match_start = next_match_start;

                    //         if next_match_start > match_start {
                    //             continue;
                    //         }
                    //     }
                    // }
                    assert!(last_match <= match_start);

                    symbols.push(Symbol::LiteralRun {
                        start: last_match as u32,
                        end: match_start as u32,
                    });

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

                    length = 0;
                    continue 'outer;
                }

                // matches.insert(current >> 8, ip + 1);
                // matches.insert(current >> 16, ip + 2);

                // If we haven't found a match in a while, start skipping ahead by emitting multiple
                // literals at once.
                ip += 1 + ((ip - last_match) >> skip_ahead_shift);
            }
            if data.len() < ip + 8 {
                symbols.push(Symbol::LiteralRun {
                    start: last_match as u32,
                    end: data.len() as u32,
                });
                ip = data.len();
            }

            Self::write_symbols(writer, data, &symbols, ip == data.len())?;
        }

        Ok(())
    }
}

const WINDOW_SIZE: usize = 32768;

struct BitWriter<W: Write> {
    buffer: u64,
    nbits: u8,
    writer: W,
}
impl<W: Write> BitWriter<W> {
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
}

/// Compressor that produces fdeflate compressed streams.
pub struct Compressor<W: Write> {
    checksum: Adler32,
    pending: Vec<u8>,
    bit_writer: BitWriter<W>,
    inner: CompressorInner,
}
impl<W: Write> Compressor<W> {
    /// Create a new Compressor.
    pub fn new(mut writer: W) -> io::Result<Self> {
        writer.write_all(&[0x78, 0x01])?; // zlib header

        let inner = match 1 {
            0 => CompressorInner::Stored,
            _ => CompressorInner::Fast {
                min_match: 4,
                skip_ahead_shift: 3,
                search_depth: 80,
                nice_length: 16,
            },
        };

        Ok(Self {
            checksum: Adler32::new(),
            bit_writer: BitWriter {
                buffer: 0,
                nbits: 0,
                writer,
            },
            pending: Vec::new(),
            inner,
        })
    }

    /// Write data to the compressor.
    pub fn write_data(&mut self, mut data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);

        // Feed the data into the compressor in 1 MB chunks.
        const CHUNK_SIZE: usize = 1 << 20;
        while !data.is_empty() {
            if self.pending.is_empty() && data.len() >= CHUNK_SIZE {
                self.inner
                    .compress_data(&mut self.bit_writer, &data[..CHUNK_SIZE], false)?;
                data = &data[CHUNK_SIZE..];
                continue;
            }

            if data.len() + self.pending.len() < CHUNK_SIZE {
                self.pending.extend_from_slice(data);
                break;
            }

            let remaining = CHUNK_SIZE - self.pending.len();
            self.pending.extend_from_slice(&data[..remaining]);
            self.inner
                .compress_data(&mut self.bit_writer, &self.pending, false)?;
            self.pending.clear();
        }

        Ok(())
    }

    /// Write the remainder of the stream and return the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.inner
            .compress_data(&mut self.bit_writer, &self.pending, true)?;

        // Write end of block
        self.bit_writer.flush()?;

        // Write Adler32 checksum
        let checksum: u32 = self.checksum.finish();
        self.bit_writer
            .writer
            .write_all(checksum.to_be_bytes().as_ref())
            .unwrap();
        Ok(self.bit_writer.writer)
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
    let mut compressed = compressor.finish().unwrap();

    if compressed.len() > StoredOnlyCompressor::<io::Cursor<&[u8]>>::compressed_size(input.len()) {
        compressed.clear();
        let mut compressor = StoredOnlyCompressor::new(io::Cursor::new(compressed)).unwrap();
        compressor.write_data(input).unwrap();
        compressor.finish().unwrap().into_inner()
    } else {
        compressed
    }
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
