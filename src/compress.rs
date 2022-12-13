use simd_adler32::Adler32;
use std::{
    convert::TryInto,
    io::{self, Write},
};

use crate::tables::{BITMASKS, HUFFMAN_CODES, HUFFMAN_LENGTHS, LEN_EXTRA, LEN_SYM};

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
            let sym = LEN_SYM[run as usize - 3] as usize;
            self.write_bits(HUFFMAN_CODES[sym] as u64, HUFFMAN_LENGTHS[sym])?;

            let len_extra = LEN_EXTRA[run as usize - 3];
            let extra = ((run - 3) & BITMASKS[len_extra as usize]) as u64;
            self.write_bits(extra, len_extra + 1)?;
        } else {
            debug_assert_eq!(HUFFMAN_CODES[0], 0);
            self.write_bits(0, run as u8 * HUFFMAN_LENGTHS[0])?;
        }

        Ok(())
    }

    pub fn new(writer: W) -> Self {
        Self {
            checksum: Adler32::new(),
            buffer: 0,
            nbits: 0,
            writer,
        }
    }

    pub fn write_headers(&mut self) -> io::Result<()> {
        self.write_bits(0x0178, 16)?; // zlib header

        self.write_bits(0b1, 1)?; // BFINAL
        self.write_bits(0b10, 2)?; // Dynamic Huffman block

        self.write_bits((HUFFMAN_LENGTHS.len() - 257) as u64, 5)?; // # of length / literal codes
        self.write_bits(0, 5)?; // 1 distance code
        self.write_bits(15, 4)?; // 16 code length codes

        // Write code lengths for code length alphabet
        for _ in 0..3 {
            self.write_bits(0, 3)?;
        }
        for _ in 0..16 {
            self.write_bits(4, 3)?;
        }

        // Write code lengths for length/literal alphabet
        for &len in &HUFFMAN_LENGTHS {
            self.write_bits((len.reverse_bits() >> 4) as u64, 4)?;
        }

        // Write code lengths for distance alphabet
        for _ in 0..1 {
            self.write_bits(0b1000, 4)?;
        }

        Ok(())
    }

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

pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new(Vec::with_capacity(input.len() / 4));
    compressor.write_headers().unwrap();
    compressor.write_data(input).unwrap();
    compressor.finish().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[bench]
    fn bench_compute_code_lengths(b: &mut test::Bencher) {
        const N: usize = 48;
        let mut rng = rand::thread_rng();
        let mut freqs = vec![0; N];
        for i in 0..freqs.len() {
            freqs[i] = rng.gen_range::<u64, _>(1..1000);
        }

        b.iter(|| {
            let mut lengths = vec![0; N];
            crate::compute_code_lengths(&freqs, &[1; N], &[8; N], &mut lengths);
        });
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

    #[bench]
    fn bench_uniform_random(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = rng.gen();
        }
        b.bytes = data.len() as u64;
        b.iter(|| compress_to_vec(&data));
    }

    #[bench]
    fn bench_low(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = (rng.gen_range::<u8, _>(0..16) * 2).wrapping_sub(16);
        }
        b.bytes = data.len() as u64;
        b.iter(|| compress_to_vec(&data));
    }

    #[bench]
    fn bench_mixture(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            if rng.gen_range(0..200) == 1 {
                *byte = rng.gen();
            } else {
                *byte = rng.gen_range::<u8, _>(0..32).wrapping_sub(16);
            }
        }
        b.bytes = data.len() as u64;
        b.iter(|| compress_to_vec(&data));
    }

    #[bench]
    fn bench_distribution(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = match rng.gen_range(0..100) {
                0 => rng.gen(),
                1..=2 => rng.gen_range::<u8, _>(0..32).wrapping_sub(16),
                11..=50 => rng.gen_range::<u8, _>(0..16).wrapping_sub(8),
                51..=80 => rng.gen_range::<u8, _>(0..8).wrapping_sub(4),
                _ => 0,
            }
        }
        b.bytes = data.len() as u64;
        b.iter(|| compress_to_vec(&data));
    }
}
