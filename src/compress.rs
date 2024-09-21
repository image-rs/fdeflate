use simd_adler32::Adler32;
use std::{collections::BinaryHeap, io::{self, Seek, SeekFrom, Write}};

use crate::tables::{
    BITMASKS, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, HUFFMAN_CODES, HUFFMAN_LENGTHS,
    LENGTH_TO_LEN_EXTRA, LENGTH_TO_SYMBOL,
};

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

    fn write_run(&mut self, mut run: u32) -> io::Result<()> {
        self.write_bits(HUFFMAN_CODES[0] as u64, HUFFMAN_LENGTHS[0])?;
        run -= 1;

        while run >= 258 {
            self.write_bits(HUFFMAN_CODES[285] as u64, HUFFMAN_LENGTHS[285] + 1)?;
            run -= 258;
        }

        if run > 4 {
            let sym = LENGTH_TO_SYMBOL[run as usize - 3] as usize;
            self.write_bits(HUFFMAN_CODES[sym] as u64, HUFFMAN_LENGTHS[sym])?;

            let len_extra = LENGTH_TO_LEN_EXTRA[run as usize - 3];
            let extra = ((run - 3) & BITMASKS[len_extra as usize]) as u64;
            self.write_bits(extra, len_extra + 1)?;
        } else {
            debug_assert_eq!(HUFFMAN_CODES[0], 0);
            self.write_bits(0, run as u8 * HUFFMAN_LENGTHS[0])?;
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

    fn match_length(data: &[u8], a: usize, b: usize) -> usize {
        if b - a > 32768 {
            return 0;
        }

        let mut length = 0;
        while length < 258 && b + length < data.len() && data[a + length] == data[b + length] {
            length += 1;
        }
        length
    }

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);

        const TABLE_SIZE: usize = 32768;
        let mut matches = vec![[0; 2]; TABLE_SIZE];

        let mut i = 0;
        while i + 8 < data.len() {

            if data[i] == 0 {
                let mut run_length = 1;
                while run_length < 258 && i + run_length < data.len() && data[i + run_length] == 0 {
                    run_length += 1;
                }
                if run_length >= 4 {
                    innumerable::event!("run", run_length.min(10) as u64);
                    self.write_run(run_length as u32)?;
                    i += run_length;
                    continue;
                }
            }

            let current = u64::from_le_bytes(data[i..][..8].try_into().unwrap());
            let current_hash = (0x27220a95u64.wrapping_mul(current & 0xffffffff)) as usize % TABLE_SIZE;
            let current_hash2 = (0x330698ec124c97f2u64.wrapping_mul(current)) as usize % TABLE_SIZE;

            let [a, b] = matches[current_hash];
            let a_len = Self::match_length(data, a, i);
            let b_len = Self::match_length(data, b, i);

            let [a2, b2] = matches[current_hash2];
            let a_len2 = Self::match_length(data, a2, i);
            let b_len2 = Self::match_length(data, b2, i);

            if a < b {
                matches[current_hash] = [i, b];
            } else {
                matches[current_hash] = [a, i];
            }
            if a2 < b2 {
                matches[current_hash2] = [i, b2];
            } else {
                matches[current_hash2] = [a2, i];
            }

            let (mut length, mut prev_i) = (a_len, a);
            if b_len > length || b_len == length && b > prev_i {
                length = b_len;
                prev_i = b;
            }
            if a_len2 > length || a_len2 == length && a2 > prev_i {
                length = a_len2;
                prev_i = a2;
            }
            if b_len2 > length || b_len2 == length && b2 > prev_i {
                length = b_len2;
                prev_i = b2;
            }

            if length >= 3 {
                let next = u64::from_le_bytes(data[i + 1..][..8].try_into().unwrap());
                let next_hash = (0x27220a95u64.wrapping_mul(next & 0xffffffff)) as usize % TABLE_SIZE;
                let next_hash2 = (0x330698ec124c97f2u64.wrapping_mul(next)) as usize % TABLE_SIZE;

                let [next_a, next_b] = matches[next_hash];
                let next_a_len = Self::match_length(data, next_a, i + 1);
                let next_b_len = Self::match_length(data, next_b, i + 1);

                let [next_a2, next_b2] = matches[next_hash2];
                let next_a_len2 = Self::match_length(data, next_a2, i + 1);
                let next_b_len2 = Self::match_length(data, next_b2, i + 1);
                let next_length = next_a_len.max(next_b_len).max(next_a_len2).max(next_b_len2);

                if length >= next_length && next != 0 {
                    let sym = LENGTH_TO_SYMBOL[length - 3] as usize;
                    let len_bits = HUFFMAN_LENGTHS[sym];
                    let len_extra = LENGTH_TO_LEN_EXTRA[length - 3];

                    let dist = (i - prev_i) as u16;
                    let mut dist_sym = 29;
                    while dist_sym > 0 && dist < DIST_SYM_TO_DIST_BASE[dist_sym as usize - 1] {
                        dist_sym -= 1;
                    }
                    let dist_bits = 6;
                    let dist_extra = DIST_SYM_TO_DIST_EXTRA[dist_sym as usize];

                    let backref_cost =
                        len_bits as u32 + len_extra as u32 + dist_bits as u32 + dist_extra as u32;
                    assert!(backref_cost < 256);
                    let backref_cost = backref_cost as u8;

                    let mut literal_cost = 0;
                    for j in i..i + length {
                        literal_cost += HUFFMAN_LENGTHS[data[j] as usize] as u32;
                    }

                    if (backref_cost as u32) < literal_cost {
                        innumerable::event!("backref", length.min(10) as u64);

                        self.write_bits(0, backref_cost)?;
                        i += length;
                        continue;
                    }
                }
            }

            innumerable::event!("literal");
            self.write_bits(
                HUFFMAN_CODES[data[i] as usize] as u64,
                HUFFMAN_LENGTHS[data[i] as usize],
            )?;
            i += 1;
        }

        for i in i..data.len() {
            innumerable::event!("literal");
            self.write_bits(
                HUFFMAN_CODES[data[i] as usize] as u64,
                HUFFMAN_LENGTHS[data[i] as usize],
            )?;
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
