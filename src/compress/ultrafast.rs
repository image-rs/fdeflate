use simd_adler32::Adler32;
use std::io::{self, Write};

use crate::tables::{
    BITMASKS, HUFFMAN_CODES, HUFFMAN_LENGTHS, LENGTH_TO_LEN_EXTRA, LENGTH_TO_SYMBOL,
};

/// Very fast zlib compressor that trades compression ratio for speed.
///
/// This compressor is designed to be fast and efficient for filtered PNG data pixel data, where it
/// is expected that there will be many long runs of zeros, and the rest of the data is mostly small
/// differences from the previous pixel. On data data that does not match this pattern, it may
/// produce output that is *larger* than the input.
pub struct UltraFastCompressor<W: Write> {
    checksum: Adler32,
    buffer: u64,
    nbits: u8,
    writer: W,
}
impl<W: Write> UltraFastCompressor<W> {
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

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);

        let mut run = 0;
        let mut chunks = data.chunks_exact(8);
        for chunk in &mut chunks {
            let ichunk = u64::from_le_bytes(chunk.try_into().unwrap());

            if ichunk == 0 {
                run += 8;
                continue;
            } else if run > 0 {
                let run_extra = ichunk.trailing_zeros() / 8;
                self.write_run(run + run_extra)?;
                run = 0;

                if run_extra > 0 {
                    run = ichunk.leading_zeros() / 8;
                    for &b in &chunk[run_extra as usize..8 - run as usize] {
                        self.write_bits(
                            HUFFMAN_CODES[b as usize] as u64,
                            HUFFMAN_LENGTHS[b as usize],
                        )?;
                    }
                    continue;
                }
            }

            let run_start = ichunk.leading_zeros() / 8;
            if run_start > 0 {
                for &b in &chunk[..8 - run_start as usize] {
                    self.write_bits(
                        HUFFMAN_CODES[b as usize] as u64,
                        HUFFMAN_LENGTHS[b as usize],
                    )?;
                }
                run = run_start;
                continue;
            }

            let n0 = HUFFMAN_LENGTHS[chunk[0] as usize];
            let n1 = HUFFMAN_LENGTHS[chunk[1] as usize];
            let n2 = HUFFMAN_LENGTHS[chunk[2] as usize];
            let n3 = HUFFMAN_LENGTHS[chunk[3] as usize];
            let bits = HUFFMAN_CODES[chunk[0] as usize] as u64
                | ((HUFFMAN_CODES[chunk[1] as usize] as u64) << n0)
                | ((HUFFMAN_CODES[chunk[2] as usize] as u64) << (n0 + n1))
                | ((HUFFMAN_CODES[chunk[3] as usize] as u64) << (n0 + n1 + n2));
            self.write_bits(bits, n0 + n1 + n2 + n3)?;

            let n4 = HUFFMAN_LENGTHS[chunk[4] as usize];
            let n5 = HUFFMAN_LENGTHS[chunk[5] as usize];
            let n6 = HUFFMAN_LENGTHS[chunk[6] as usize];
            let n7 = HUFFMAN_LENGTHS[chunk[7] as usize];
            let bits2 = HUFFMAN_CODES[chunk[4] as usize] as u64
                | ((HUFFMAN_CODES[chunk[5] as usize] as u64) << n4)
                | ((HUFFMAN_CODES[chunk[6] as usize] as u64) << (n4 + n5))
                | ((HUFFMAN_CODES[chunk[7] as usize] as u64) << (n4 + n5 + n6));
            self.write_bits(bits2, n4 + n5 + n6 + n7)?;
        }

        if run > 0 {
            self.write_run(run)?;
        }

        for &b in chunks.remainder() {
            self.write_bits(
                HUFFMAN_CODES[b as usize] as u64,
                HUFFMAN_LENGTHS[b as usize],
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
        let mut compressor = UltraFastCompressor::new(Vec::with_capacity(input.len() / 4)).unwrap();
        compressor.write_data(input).unwrap();
        compressor.finish().unwrap()
    }

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
