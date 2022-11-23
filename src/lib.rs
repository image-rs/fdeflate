#![feature(stdsimd)]
#![feature(test)]
extern crate test;

use std::io::{self, Cursor, Write};

use simd_adler32::Adler32;

pub mod decode;

#[derive(Copy, Clone, Debug)]
pub enum Flush {
    Normal,
    Finish,
}

pub enum Status {
    Done,
    NeedMoreOutputSpace,
}

/// Length code for length values.
pub const LEN_SYM: [u16; 256] = [
    257, 258, 259, 260, 261, 262, 263, 264, 265, 265, 266, 266, 267, 267, 268, 268, 269, 269, 269,
    269, 270, 270, 270, 270, 271, 271, 271, 271, 272, 272, 272, 272, 273, 273, 273, 273, 273, 273,
    273, 273, 274, 274, 274, 274, 274, 274, 274, 274, 275, 275, 275, 275, 275, 275, 275, 275, 276,
    276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277,
    277, 277, 277, 277, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278,
    278, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 280, 280,
    280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 281, 281, 281, 281, 281,
    281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281,
    281, 281, 281, 281, 281, 281, 281, 281, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282,
    282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282,
    282, 282, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283,
    283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 284, 284, 284, 284,
    284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284,
    284, 284, 284, 284, 284, 284, 284, 284, 285,
];

/// Number of extra bits for length values.
pub const LEN_EXTRA: [u8; 256] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0,
];

const BITMASKS: [u32; 17] = [
    0x0000, 0x0001, 0x0003, 0x0007, 0x000F, 0x001F, 0x003F, 0x007F, 0x00FF, 0x01FF, 0x03FF, 0x07FF,
    0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF,
];

const HUFFMAN_LENGTHS: [u8; 286] = [
    2, 3, 4, 5, 5, 6, 6, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10,
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 11, 11, 11, 11,
    11, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 7,
    7, 7, 6, 6, 6, 5, 4, 3, 12, 8, 8, 9, 9, 11, 10, 11, 11, 10, 11, 11, 11, 11, 11, 11, 12, 11, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 9,
];

// Dynamic programming huffman tree algorithm from fpnge.
pub fn compute_code_lengths(
    freqs: &[u64],
    min_limit: &[u8],
    max_limit: &[u8],
    calculated_nbits: &mut [u8],
) {
    debug_assert_eq!(freqs.len(), min_limit.len());
    debug_assert_eq!(freqs.len(), max_limit.len());
    debug_assert_eq!(freqs.len(), calculated_nbits.len());
    let len = freqs.len();

    for i in 0..len {
        debug_assert!(min_limit[i] >= 1);
        debug_assert!(min_limit[i] <= max_limit[i]);
    }

    let precision = *max_limit.iter().max().unwrap();
    let num_patterns = 1 << precision;

    let mut dynp = vec![u64::MAX; (num_patterns + 1) * (len + 1)];
    let index = |sym: usize, off: usize| sym * (num_patterns + 1) + off;

    dynp[index(0, 0)] = 0;
    for sym in 0..len {
        for bits in min_limit[sym]..=max_limit[sym] {
            let off_delta = 1 << (precision - bits);
            for off in 0..=num_patterns.saturating_sub(off_delta) {
                dynp[index(sym + 1, off + off_delta)] = dynp[index(sym, off)]
                    .saturating_add(freqs[sym] * u64::from(bits))
                    .min(dynp[index(sym + 1, off + off_delta)]);
            }
        }
    }

    let mut sym = len;
    let mut off = num_patterns;

    while sym > 0 {
        sym -= 1;
        assert!(off > 0);

        for bits in min_limit[sym]..=max_limit[sym] {
            let off_delta = 1 << (precision - bits);
            if off_delta <= off
                && dynp[index(sym + 1, off)]
                    == dynp[index(sym, off - off_delta)]
                        .saturating_add(freqs[sym] * u64::from(bits))
            {
                off -= off_delta;
                calculated_nbits[sym] = bits;
                break;
            }
        }
    }

    for i in 0..len {
        debug_assert!(calculated_nbits[i] >= min_limit[i]);
        debug_assert!(calculated_nbits[i] <= max_limit[i]);
    }
}

fn compute_codes<const NSYMS: usize>(lengths: &[u8; NSYMS]) -> [u16; NSYMS] {
    let mut codes = [0u16; NSYMS];

    let max_length = *lengths.iter().max().unwrap();

    let mut code = 0u16;
    for len in 1..=max_length {
        for (i, &length) in lengths.iter().enumerate() {
            if length == len {
                codes[i] = code.reverse_bits() >> (16 - len);
                code += 1;
            }
        }
        code <<= 1;
    }

    assert_eq!(code, 2 << max_length);

    codes
}

pub struct Compressor<W: Write> {
    buffer: u64,
    nbits: u8,
    bytes_written: usize,

    codes: [u16; 286],
    lengths: [u8; 286],

    writer: W,
    checksum: Adler32,
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
            self.bytes_written += 8;
        }
        debug_assert!(self.nbits < 64);
        Ok(())
    }

    fn flush_to_byte_boundary(&mut self) -> io::Result<()> {
        if self.nbits % 8 != 0 {
            self.write_bits(0, 8 - self.nbits % 8)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.flush_to_byte_boundary()?;
        if self.nbits > 0 {
            self.writer
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])
                .unwrap();
            self.bytes_written += self.nbits as usize / 8;
            self.buffer = 0;
            self.nbits = 0;
        }
        Ok(())
    }

    fn write_run(&mut self, mut run: u32) -> io::Result<()> {
        self.write_bits(self.codes[0] as u64, self.lengths[0])?;
        run -= 1;

        while run >= 258 {
            self.write_bits(self.codes[285] as u64, self.lengths[285] + 1)?;
            run -= 258;
        }

        if run >= 3 {
            let sym = LEN_SYM[run as usize - 3] as usize;
            self.write_bits(self.codes[sym] as u64, self.lengths[sym])?;

            let len_extra = LEN_EXTRA[run as usize - 3];
            let extra = ((run - 3) & BITMASKS[len_extra as usize]) as u64;
            self.write_bits(extra, len_extra + 1)?;
            run = 0;
        }

        for _ in 0..run {
            self.write_bits(self.codes[0] as u64, self.lengths[0])?;
        }

        Ok(())
    }

    fn compress_simple(&mut self, data: &[u8]) -> io::Result<()> {
        self.checksum.write(data);

        const CHUNK_SIZE: usize = 64;

        let mut run = 0;
        let mut chunks = data.chunks_exact(CHUNK_SIZE);
        for chunk in &mut chunks {
            let mut mask = 0u64;
            for s in chunk.chunks_exact(8) {
                mask |= u64::from_ne_bytes(s.try_into().unwrap());
            }
            if mask == 0 {
                run += CHUNK_SIZE as u32;
                continue;
            }

            for chunk in chunk.chunks_exact(4) {
                if chunk == [0; 4] {
                    run += 4;
                    continue;
                } else if run > 0 {
                    self.write_run(run)?;
                    run = 0;
                }

                let n0 = self.lengths[chunk[0] as usize];
                let n1 = self.lengths[chunk[1] as usize];
                let n2 = self.lengths[chunk[2] as usize];
                let n3 = self.lengths[chunk[3] as usize];

                let bits = self.codes[chunk[0] as usize] as u64
                    | ((self.codes[chunk[1] as usize] as u64) << n0)
                    | ((self.codes[chunk[2] as usize] as u64) << (n0 + n1))
                    | ((self.codes[chunk[3] as usize] as u64) << (n0 + n1 + n2));

                let nbits = n0 + n1 + n2 + n3;
                self.write_bits(bits, nbits)?;
            }
        }

        if run > 0 {
            self.write_run(run)?;
        }

        for &b in chunks.remainder() {
            self.write_bits(self.codes[b as usize] as u64, self.lengths[b as usize])?;
        }

        Ok(())
    }

    // fn compress_simple2(
    //     &mut self,
    //     data: &[u8],
    //     lengths: &[u8],
    //     codes: &[u16],
    //     out_buf: &mut &mut [u8],
    // ) {
    //     let mut run = 0;
    //     let mut chunks = data.chunks_exact(24);
    //     for chunk in &mut chunks {
    //         if chunk == [0; 24] {
    //             run += 24;
    //         } else if chunk.chunks_exact(4).any(|c| c == [0; 4]) {
    //             for chunk in chunk.chunks_exact(4) {
    //                 if chunk == [0; 4] {
    //                     run += 4;
    //                     continue;
    //                 } else if run > 0 {
    //                     self.write_run(run, lengths, codes, out_buf);
    //                     run = 0;
    //                 }

    //                 let n0 = lengths[chunk[0] as usize];
    //                 let n1 = lengths[chunk[1] as usize];
    //                 let n2 = lengths[chunk[2] as usize];
    //                 let n3 = lengths[chunk[3] as usize];

    //                 let bits = codes[chunk[0] as usize] as u64
    //                     | ((codes[chunk[1] as usize] as u64) << n0)
    //                     | ((codes[chunk[2] as usize] as u64) << (n0 + n1))
    //                     | ((codes[chunk[3] as usize] as u64) << (n0 + n1 + n2));

    //                 let nbits = n0 + n1 + n2 + n3;
    //                 self.write_bits(bits, nbits, out_buf);
    //             }
    //         } else {
    //             for chunk in chunk.chunks_exact(4) {
    //                 let n0 = lengths[chunk[0] as usize];
    //                 let n1 = lengths[chunk[1] as usize];
    //                 let n2 = lengths[chunk[2] as usize];
    //                 let n3 = lengths[chunk[3] as usize];

    //                 let bits = codes[chunk[0] as usize] as u64
    //                     | ((codes[chunk[1] as usize] as u64) << n0)
    //                     | ((codes[chunk[2] as usize] as u64) << (n0 + n1))
    //                     | ((codes[chunk[3] as usize] as u64) << (n0 + n1 + n2));

    //                 let nbits = n0 + n1 + n2 + n3;
    //                 self.write_bits(bits, nbits, out_buf);
    //             }
    //         }
    //     }

    //     if run > 0 {
    //         self.write_run(run, lengths, codes, out_buf);
    //     }

    //     for &b in chunks.remainder() {
    //         self.write_bits(codes[b as usize] as u64, lengths[b as usize], out_buf);
    //     }
    // }

    // fn compress_simple3(
    //     &mut self,
    //     data: &[u8],
    //     lengths: &[u8],
    //     codes: &[u16],
    //     out_buf: &mut &mut [u8],
    // ) {
    //     let mut run = 0;
    //     let mut chunks = data.chunks_exact(8);
    //     for chunk in &mut chunks {
    //         if chunk == [0; 8] {
    //             run += 8;
    //             continue;
    //         } else if run > 0 {
    //             self.write_run(run, lengths, codes, out_buf);
    //             run = 0;
    //         }

    //         let n0 = lengths[chunk[0] as usize];
    //         let n1 = lengths[chunk[1] as usize];
    //         let n2 = lengths[chunk[2] as usize];
    //         let n3 = lengths[chunk[3] as usize];

    //         let bits = codes[chunk[0] as usize] as u64
    //             | ((codes[chunk[1] as usize] as u64) << n0)
    //             | ((codes[chunk[2] as usize] as u64) << (n0 + n1))
    //             | ((codes[chunk[3] as usize] as u64) << (n0 + n1 + n2));

    //         let nbits = n0 + n1 + n2 + n3;
    //         self.write_bits(bits, nbits, out_buf);

    //         let n0 = lengths[chunk[4] as usize];
    //         let n1 = lengths[chunk[5] as usize];
    //         let n2 = lengths[chunk[6] as usize];
    //         let n3 = lengths[chunk[7] as usize];

    //         let bits = codes[chunk[4] as usize] as u64
    //             | ((codes[chunk[5] as usize] as u64) << n0)
    //             | ((codes[chunk[6] as usize] as u64) << (n0 + n1))
    //             | ((codes[chunk[7] as usize] as u64) << (n0 + n1 + n2));

    //         let nbits = n0 + n1 + n2 + n3;
    //         self.write_bits(bits, nbits, out_buf);
    //     }

    //     if run > 0 {
    //         self.write_run(run, lengths, codes, out_buf);
    //     }

    //     for &b in chunks.remainder() {
    //         self.write_bits(codes[b as usize] as u64, lengths[b as usize], out_buf);
    //     }
    // }

    pub fn new(writer: W) -> Self {
        let codes = compute_codes(&HUFFMAN_LENGTHS);

        Self {
            buffer: 0,
            nbits: 0,
            bytes_written: 0,
            codes,
            lengths: HUFFMAN_LENGTHS,
            writer,
            checksum: Adler32::new(),
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

    pub fn write_data(&mut self, in_buf: &[u8]) -> io::Result<()> {
        self.compress_simple(in_buf)
    }

    pub fn finish(mut self) -> io::Result<W> {
        // Write end of block
        self.write_bits(self.codes[256] as u64, HUFFMAN_LENGTHS[256])?;
        self.flush()?;

        // Write Adler32 checksum
        let checksum: u32 = self.checksum.finish();
        self.writer
            .write_all(checksum.to_be_bytes().as_ref())
            .unwrap();
        self.bytes_written += 4;
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
            compute_code_lengths(&freqs, &[1; N], &[8; N], &mut lengths);
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
            //*byte = (rng.gen_range::<u8, _>(0..16) * 2).wrapping_sub(16);
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
