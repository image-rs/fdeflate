#![feature(stdsimd)]
#![feature(test)]
extern crate test;

use std::io::Write;

#[derive(Copy, Clone, Debug)]
pub enum Flush {
    Normal,
    Finish,
}

pub enum Status {
    Done,
    NeedMoreOutputSpace,
}

const BLOCK_SIZE: usize = 24;

/// Length code for length values.
#[rustfmt::skip]
pub const LEN_SYM: [u16; 256] = [
    257, 258, 259, 260, 261, 262, 263, 264, 265, 265, 266, 266, 267, 267, 268, 268,
    269, 269, 269, 269, 270, 270, 270, 270, 271, 271, 271, 271, 272, 272, 272, 272,
    273, 273, 273, 273, 273, 273, 273, 273, 274, 274, 274, 274, 274, 274, 274, 274,
    275, 275, 275, 275, 275, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 276,
    277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277,
    278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278, 278,
    279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279, 279,
    280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280, 280,
    281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281,
    281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281, 281,
    282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282,
    282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282, 282,
    283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283,
    283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283, 283,
    284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284,
    284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 284, 285
];

/// Number of extra bits for length values.
#[rustfmt::skip]
pub const LEN_EXTRA: [u8; 256] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0
];

#[rustfmt::skip]
const BITMASKS: [u32; 17] = [
    0x0000, 0x0001, 0x0003, 0x0007, 0x000F, 0x001F, 0x003F, 0x007F, 0x00FF,
    0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF
];

// const DYN3_HUFF_LENGTHS: [u8; 288] = [
//     2, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
//     10, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 12, 11, 12, 12, 11, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
//     12, 12, 12, 12, 12, 12, 11, 12, 12, 11, 12, 11, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
//     10, 9, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5,
//     5, 4, 12, 6, 0, 0, 8, 0, 0, 8, 0, 10, 0, 10, 12, 10, 12, 12, 12, 12, 12, 12, 12, 12, 9, 12, 12,
//     12, 10, 12, 6, 12, 0, 0,
// ];

const DYN3_HUFF_LENGTHS: [u8; 288] = [
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
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 9, 0, 0,
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

// fn build_literal_length_tree(buf: &[u8]) -> Vec<u8> {
//     let mut counts = vec![1u64; 47];
//     counts[46] = 1;
//     let sample_length = 2048.min(buf.len());
//     for &b in &buf[..sample_length] {
//         match b {
//             0..=15 => counts[b as usize] += 1,
//             240..=255 => counts[b as usize - 240 + 16] += 1,
//             b => counts[((b - 16) >> 4) as usize + 32] += 1,
//         }
//     }

//     let mut min_limit = vec![1u8; 47];
//     let max_limit = vec![8u8; 47];
//     // for i in 32..46 {
//     //     min_limit[i] = 8;
//     // }

//     let mut lengths = vec![0u8; 47];
//     compute_code_lengths(&counts, &min_limit, &max_limit, &mut lengths);

//     let mut output_lengths = Vec::with_capacity(257);
//     output_lengths.extend_from_slice(&lengths[..16]);
//     for &len in &lengths[32..46] {
//         output_lengths.extend_from_slice(&[len + 4; 16]);
//     }
//     output_lengths.extend_from_slice(&lengths[16..32]);
//     output_lengths.push(lengths[46]);
//     assert_eq!(output_lengths.len(), 257);

//     // for i in 16..240 {
//     //     assert_eq!(output_lengths[i], 12);
//     // }

//     output_lengths
// }

// fn build_literal_length_tree(buf: &[u8]) -> Vec<u8> {
//     let mut counts = vec![1u64; 257];
//     let sample_length = 2048000.min(buf.len());
//     for &b in &buf[..sample_length] {
//         counts[b as usize] += 1;
//     }
//     let min_limit = vec![1u8; 257];
//     let mut max_limit = vec![12u8; 257];
//     for i in 0..16 {
//         max_limit[i] = 8;
//         max_limit[240 + i] = 8;
//     }
//     let mut lengths = vec![0u8; 257];
//     compute_code_lengths(&counts, &min_limit, &max_limit, &mut lengths);
//     lengths
// }

fn build_literal_length_tree(_buf: &[u8]) -> Vec<u8> {
    DYN3_HUFF_LENGTHS.to_vec()
}

fn compute_codes(lengths: &[u8]) -> Vec<u16> {
    let mut codes = vec![0u16; lengths.len()];

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

// #[repr(align(32))]
// struct AlignedBlock([u8; BLOCK_SIZE]);

// #[repr(align(32))]
// struct AlignedBlockU64([u64; BLOCK_SIZE / 8]);

pub struct Compressor {
    buffer: u64,
    nbits: u8,
    bytes_written: usize,
}
impl Compressor {
    pub fn new() -> Compressor {
        Compressor {
            buffer: 0,
            nbits: 0,
            bytes_written: 0,
        }
    }

    fn write_bits(&mut self, bits: u64, nbits: u8, output: &mut &mut [u8]) {
        debug_assert!(nbits <= 64);

        self.buffer |= bits << self.nbits;
        self.nbits += nbits;

        if self.nbits >= 64 {
            output.write_all(&self.buffer.to_le_bytes()).unwrap();
            self.nbits -= 64;
            self.buffer = bits.checked_shr((nbits - self.nbits) as u32).unwrap_or(0);
            self.bytes_written += 8;
        }
        debug_assert!(self.nbits < 64);
    }

    fn flush_to_byte_boundary(&mut self, output: &mut &mut [u8]) {
        if self.nbits % 8 != 0 {
            self.write_bits(0, 8 - self.nbits % 8, output);
        }
    }

    fn flush(&mut self, output: &mut &mut [u8]) {
        self.flush_to_byte_boundary(output);
        if self.nbits > 0 {
            output
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])
                .unwrap();
            self.bytes_written += self.nbits as usize / 8;
            self.buffer = 0;
            self.nbits = 0;
        }
    }

    fn compress_aligned(
        &mut self,
        data: &[u8],
        lengths: &[u8],
        codes: &[u16],
        out_buf: &mut &mut [u8],
    ) {
        // let mut codes_low = AlignedBlock([0; BLOCK_SIZE]);
        // let mut codes_high = AlignedBlock([0; BLOCK_SIZE]);
        // let mut lengths_low = AlignedBlock([0; BLOCK_SIZE]);
        // let mut lengths_high = AlignedBlock([0; BLOCK_SIZE]);
        // for i in 0..16 {
        //     assert!(lengths[i] <= 8);
        //     codes_low.0[i] = codes[i] as u8;
        //     lengths_low.0[i] = lengths[i];
        //     if BLOCK_SIZE == 32 {
        //         codes_low.0[i + 16] = codes[i] as u8;
        //         lengths_low.0[i + 16] = lengths[i];
        //     }
        // }
        // for i in 0..16 {
        //     assert!(lengths[i] <= 8);
        //     codes_high.0[i] = codes[i + 240] as u8;
        //     lengths_high.0[i] = lengths[i + 240];
        //     if BLOCK_SIZE == 32 {
        //         codes_high.0[i + 16] = codes[i + 240] as u8;
        //         lengths_high.0[i + 16] = lengths[i + 240];
        //     }
        // }

        let mut short_codes = [0; 32];
        let mut short_lengths = [0; 32];
        for i in 0..16 {
            // short_codes[i + 16] = codes[i] as u8;
            // short_lengths[i + 16] = lengths[i];
            // short_codes[i] = codes[i + 240] as u8;
            // short_lengths[i] = lengths[i + 240];
            short_codes[i] = codes[i] as u8;
            short_lengths[i] = lengths[i];
            short_codes[i + 16] = codes[i + 240] as u8;
            short_lengths[i + 16] = lengths[i + 240];
        }

        let all_zero_codes = codes[0] as u64
            | (codes[0] as u64) << lengths[0]
            | (codes[0] as u64) << (lengths[0] * 2)
            | (codes[0] as u64) << (lengths[0] * 3)
            | (codes[0] as u64) << (lengths[0] * 4)
            | (codes[0] as u64) << (lengths[0] * 5)
            | (codes[0] as u64) << (lengths[0] * 6)
            | (codes[0] as u64) << (lengths[0] * 7);
        let all_zero_length = 8 * lengths[0];

        let zero_repeat24_code = codes[0] as u64
            | ((codes[270] as u64) << lengths[0])
            | (0b00 << (lengths[0] + lengths[270]));
        let zero_repeat24_length = lengths[0] + lengths[270] + 2 + 1;

        /* unsafe */
        {
            // let codes_low = _mm256_load_si256(codes_low.0.as_ptr() as *const _);
            // let codes_high = _mm256_load_si256(codes_high.0.as_ptr() as *const _);
            // let lengths_low = _mm256_load_si256(lengths_low.0.as_ptr() as *const _);
            // let lengths_high = _mm256_load_si256(lengths_high.0.as_ptr() as *const _);

            // let codes_low = _mm_load_si128(codes_low.0.as_ptr() as *const _);
            // let codes_high = _mm_load_si128(codes_high.0.as_ptr() as *const _);
            // let lengths_low = _mm_load_si128(lengths_low.0.as_ptr() as *const _);
            // let lengths_high = _mm_load_si128(lengths_high.0.as_ptr() as *const _);

            let mut run = 0;
            assert_eq!(data.len() % BLOCK_SIZE, 0);
            for chunk in data.chunks_exact(BLOCK_SIZE) {
                // let d = _mm256_load_si256(chunk.as_ptr() as *const _);
                // let use_lohi = _mm256_cmpgt_epi8(
                //     _mm256_add_epi8(d, _mm256_set1_epi8(112)),
                //     _mm256_set1_epi8(95),
                // );
                // let use_lohi_mask = _mm256_movemask_epi8(use_lohi) == u32::MAX as i32;

                // let d = _mm_load_si128(chunk.as_ptr() as *const _);
                // let use_lohi =
                //     _mm_cmpgt_epi8(_mm_add_epi8(d, _mm_set1_epi8(112)), _mm_set1_epi8(95));
                // let use_lohi_mask = _mm_movemask_epi8(use_lohi) == 0xffff;

                // let use_lohi_mask = chunk.iter().all(|&b| (b as i8).wrapping_add(112) > 95);

                let mut mask = 0;
                for &b in chunk {
                    mask |= b.wrapping_add(16);
                }
                let use_lohi_mask = mask & 0xe0 == 0;

                if use_lohi_mask && false {
                    // if mask == 16 {
                    //     let mut mask = 0;
                    //     for &b in chunk {
                    //         mask |= b;
                    //     }
                    //     if mask == 0 {
                    //         for _ in 0..BLOCK_SIZE / 8 {
                    //             self.write_bits(all_zero_codes, all_zero_length, out_buf);
                    //         }
                    //         print!(".");
                    //         continue;
                    //     }
                    // }
                    // let mut bits_mem = AlignedBlockU64([0u64; BLOCK_SIZE / 8]);
                    // let mut bitmask_mem = AlignedBlockU64([0u64; BLOCK_SIZE / 8]);
                    // let mut bitcount_mem = AlignedBlockU64([0u64; BLOCK_SIZE / 8]);

                    // let data_for_lut = _mm_and_si128(d, _mm_set1_epi8(0x0F));
                    // let select_low_high = _mm_cmpgt_epi8(d, _mm_set1_epi8(16));

                    // let bits_low = _mm_shuffle_epi8(codes_low, data_for_lut);
                    // let bits_high = _mm_shuffle_epi8(codes_high, data_for_lut);
                    // let bits = _mm_blendv_epi8(bits_low, bits_high, select_low_high);
                    // _mm_store_si128(bits_mem.0.as_mut_ptr() as *mut _, bits);

                    // let nbits_low = _mm_shuffle_epi8(lengths_low, data_for_lut);
                    // let nbits_high = _mm_shuffle_epi8(lengths_high, data_for_lut);
                    // let nbits = _mm_blendv_epi8(nbits_low, nbits_high, select_low_high);

                    // let bitmask = _mm_shuffle_epi8(
                    //     _mm_set_epi32(
                    //         0xffffffffu32 as i32,
                    //         0xffffffffu32 as i32,
                    //         0x7f3f1f0f,
                    //         0x07030100,
                    //     ),
                    //     nbits,
                    // );
                    // _mm_store_si128(bitmask_mem.0.as_mut_ptr() as *mut _, bitmask);

                    // let bitcount = _mm_sad_epu8(nbits, _mm_setzero_si128());
                    // _mm_store_si128(bitcount_mem.0.as_mut_ptr() as *mut _, bitcount);

                    // for i in 0..(BLOCK_SIZE / 8) {
                    //     let b = _pext_u64(bits_mem.0[i], bitmask_mem.0[i]);
                    //     let n = bitcount_mem.0[i];
                    //     self.write_bits(b, n as u8, out_buf);
                    // }

                    // {
                    //     let bits: &mut [u8] = bytemuck::cast_slice_mut(&mut bits);
                    //     let masks: &mut [u8] = bytemuck::cast_slice_mut(&mut masks);
                    //     for (i, &b) in chunk.iter().enumerate() {
                    //         let v = b.wrapping_add(16);
                    //         bits[i] = short_codes[v as usize];
                    //         masks[i] = ((1u16 << short_lengths[v as usize]) - 1) as u8;
                    //     }
                    // }
                    // for i in 0..(BLOCK_SIZE/8) {
                    //     let b = _pext_u64(bits[i], masks[i] as u64);
                    //     let n = masks[i].count_ones();
                    //     self.write_bits(b, n as u8, out_buf);
                    // }

                    // let mut vs = [0u8; BLOCK_SIZE];
                    // _mm256_storeu_si256(vs.as_mut_ptr() as *mut _, sum);
                    // for v in vs {
                    //     self.write_bits(short_codes[v as usize] as u64, short_lengths[v as usize] as u8, out_buf);
                    // }

                    // for b in chunk {
                    //     let v = (b & 0x1f) as usize;
                    //     self.write_bits(short_codes[v] as u64, short_lengths[v] as u8, out_buf);
                    // }

                    // if chunk.iter().all(|&c| c == chunk[0]) {
                    //     continue;
                    // }

                    for c in chunk.chunks_exact(8) {
                        if c == &[0; 8] {
                            self.write_bits(all_zero_codes, all_zero_length, out_buf);
                            continue;
                        }

                        let n0 = short_lengths[c[0] as usize & 0x1f];
                        let n1 = short_lengths[c[1] as usize & 0x1f];
                        let n2 = short_lengths[c[2] as usize & 0x1f];
                        let n3 = short_lengths[c[3] as usize & 0x1f];
                        let n4 = short_lengths[c[4] as usize & 0x1f];
                        let n5 = short_lengths[c[5] as usize & 0x1f];
                        let n6 = short_lengths[c[6] as usize & 0x1f];
                        let n7 = short_lengths[c[7] as usize & 0x1f];
                        let sum4 = n0 + n1 + n2 + n3;

                        let bits = short_codes[c[0] as usize & 0x1f] as u64
                            | ((short_codes[c[1] as usize & 0x1f] as u64) << n0)
                            | ((short_codes[c[2] as usize & 0x1f] as u64) << (n0 + n1))
                            | ((short_codes[c[3] as usize & 0x1f] as u64) << (n0 + n1 + n2))
                            | ((short_codes[c[4] as usize & 0x1f] as u64) << (n0 + n1 + n2 + n3))
                            | ((short_codes[c[5] as usize & 0x1f] as u64) << (sum4 + n4))
                            | ((short_codes[c[6] as usize & 0x1f] as u64) << (sum4 + n4 + n5))
                            | ((short_codes[c[7] as usize & 0x1f] as u64) << (sum4 + n4 + n5 + n6));

                        let nbits = sum4 + n4 + n5 + n6 + n7;
                        self.write_bits(bits, nbits, out_buf);
                    }

                    continue;
                }

                if chunk == &[0; BLOCK_SIZE] {
                    self.write_bits(zero_repeat24_code, zero_repeat24_length, out_buf);
                    continue;
                }

                for c in chunk.chunks_exact(4) {
                    let n0 = lengths[c[0] as usize];
                    let n1 = lengths[c[1] as usize];
                    let n2 = lengths[c[2] as usize];
                    let n3 = lengths[c[3] as usize];

                    let bits = codes[c[0] as usize] as u64
                        | ((codes[c[1] as usize] as u64) << n0)
                        | ((codes[c[2] as usize] as u64) << (n0 + n1))
                        | ((codes[c[3] as usize] as u64) << (n0 + n1 + n2));

                    let nbits = n0 + n1 + n2 + n3;
                    self.write_bits(bits, nbits, out_buf);
                }

                // for &b in chunk {
                //     self.write_bits(codes[b as usize] as u64, lengths[b as usize], out_buf);
                // }
            }
        }
    }

    fn write_run(&mut self, mut run: u32, lengths: &[u8], codes: &[u16], out_buf: &mut &mut [u8]) {
        self.write_bits(codes[0] as u64, lengths[0], out_buf);
        run -= 1;

        while run >= 258 {
            self.write_bits(codes[285] as u64, lengths[285] + 1, out_buf);
            run -= 258;
        }

        if run >= 3 {
            let sym = LEN_SYM[run as usize - 3] as usize;
            self.write_bits(codes[sym] as u64, lengths[sym], out_buf);

            let len_extra = LEN_EXTRA[run as usize - 3];
            let extra = ((run - 3) & BITMASKS[len_extra as usize]) as u64;
            self.write_bits(extra, len_extra + 1, out_buf);
            run = 0;
        }

        for _ in 0..run {
            self.write_bits(codes[0] as u64, lengths[0], out_buf);
        }
    }

    fn compress_simple(
        &mut self,
        data: &[u8],
        lengths: &[u8],
        codes: &[u16],
        out_buf: &mut &mut [u8],
    ) {
        let mut run = 0;
        let mut chunks = data.chunks_exact(4);
        for chunk in &mut chunks {
            if chunk == &[0; 4] {
                run += 4;
                continue;
            } else if run > 0 {
                self.write_run(run, lengths, codes, out_buf);
                run = 0;
            }

            /*for chunk in chunk.chunks_exact(4)*/
            {
                let n0 = lengths[chunk[0] as usize];
                let n1 = lengths[chunk[1] as usize];
                let n2 = lengths[chunk[2] as usize];
                let n3 = lengths[chunk[3] as usize];

                let bits = codes[chunk[0] as usize] as u64
                    | ((codes[chunk[1] as usize] as u64) << n0)
                    | ((codes[chunk[2] as usize] as u64) << (n0 + n1))
                    | ((codes[chunk[3] as usize] as u64) << (n0 + n1 + n2));

                let nbits = n0 + n1 + n2 + n3;
                self.write_bits(bits, nbits, out_buf);
            }
        }

        if run > 0 {
            self.write_run(run, lengths, codes, out_buf);
        }

        for &b in chunks.remainder() {
            self.write_bits(codes[b as usize] as u64, lengths[b as usize], out_buf);
        }
    }

    /// Returns (status, bytes_in, bytes_out)
    pub fn compress(
        &mut self,
        in_buf: &[u8],
        mut out_buf: &mut [u8],
        _flush: Flush,
    ) -> (Status, usize, usize) {
        self.write_bits(0x0178, 16, &mut out_buf); // zlib header

        self.write_bits(0b1, 1, &mut out_buf); // BFINAL
        self.write_bits(0b10, 2, &mut out_buf); // Dynamic Huffman block

        let lengths = build_literal_length_tree(in_buf);
        //assert_eq!(lengths.len(), 257);

        self.write_bits((lengths.len() - 257) as u64, 5, &mut out_buf); // # of length / literal codes
        self.write_bits(0, 5, &mut out_buf); // 1 distance code
        self.write_bits(15, 4, &mut out_buf); // 19 code length codes

        // Write code lengths for code length alphabet
        for _ in 0..3 {
            self.write_bits(0, 3, &mut out_buf);
        }
        for _ in 0..16 {
            self.write_bits(4, 3, &mut out_buf);
        }

        for &len in &lengths {
            self.write_bits((len.reverse_bits() >> 4) as u64, 4, &mut out_buf);
        }

        // Write code lengths for distance alphabet
        for _ in 0..1 {
            self.write_bits(0b1000, 4, &mut out_buf); // 1 bit for distance code 0
        }

        // Write data
        let codes = compute_codes(&lengths);

        self.compress_simple(in_buf, &lengths, &codes, &mut out_buf);

        // let unaligned_prefix =
        //     ((BLOCK_SIZE - (in_buf.as_ptr() as usize % BLOCK_SIZE)) % BLOCK_SIZE).min(in_buf.len());
        // let unaligned_suffix = (in_buf.len() - unaligned_prefix) % BLOCK_SIZE;
        // let aligned_length = in_buf.len() - unaligned_prefix - unaligned_suffix;
        // assert_eq!(aligned_length % BLOCK_SIZE, 0);

        // for &byte in &in_buf[..unaligned_prefix] {
        //     self.write_bits(
        //         codes[byte as usize] as u64,
        //         lengths[byte as usize],
        //         &mut out_buf,
        //     );
        // }
        // if aligned_length > 0 {
        //     self.compress_aligned(
        //         &in_buf[unaligned_prefix..][..aligned_length],
        //         &lengths,
        //         &codes,
        //         &mut out_buf,
        //     );
        // }
        // for &byte in &in_buf[in_buf.len() - unaligned_suffix..] {
        //     self.write_bits(
        //         codes[byte as usize] as u64,
        //         lengths[byte as usize],
        //         &mut out_buf,
        //     );
        // }

        // Write end of block
        self.write_bits(codes[256] as u64, lengths[256], &mut out_buf);
        self.flush(&mut out_buf);

        // Write Adler32 checksum
        let checksum: u32 = simd_adler32::adler32(&in_buf);
        out_buf.write_all(checksum.to_be_bytes().as_ref()).unwrap();
        self.bytes_written += 4;

        (Status::Done, in_buf.len(), self.bytes_written)
    }
}

pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new();
    let mut output = vec![0; ::core::cmp::max(1024 + input.len() * 2, 2)];

    let mut in_pos = 0;
    let mut out_pos = 0;
    loop {
        let (status, bytes_in, bytes_out) =
            compressor.compress(&input[in_pos..], &mut output[out_pos..], Flush::Finish);

        out_pos += bytes_out;
        in_pos += bytes_in;

        match status {
            Status::Done => {
                output.truncate(out_pos);
                break;
            }
            Status::NeedMoreOutputSpace => {
                // We need more space, so resize the vector.
                if output.len().saturating_sub(out_pos) < 30 {
                    output.resize(output.len() * 2, 0)
                }
            }
        }
    }

    output.resize(out_pos, 0);
    output
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
