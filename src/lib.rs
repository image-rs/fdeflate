//! A fast deflate implementation.
//!
//! This crate contains an optimized implementation of the deflate algorithm tuned to compress PNG
//! images. It is compatible with standard zlib, but make a bunch of simplifying assumptions that
//! drastically improve encoding performance:
//!
//! - Exactly one block per deflate stream.
//! - No distance codes except for run length encoding of zeros.
//! - A single fixed huffman tree trained on a large corpus of PNG images.
//! - All huffman codes are 12 bits or less.
//!
//! It also contains a fast decompressor that supports arbitrary zlib streams but does especially
//! well on streams that meet the above assumptions.
//!
//! # Inspiration
//!
//! The algorithms in this crate take inspiration from multiple sources:
//! * [fpnge](https://github.com/veluca93/fpnge)
//! * [zune-inflate](https://github.com/etemesi254/zune-image/tree/main/zune-inflate)
//! * [RealTime Data Compression blog](https://fastcompression.blogspot.com/2015/10/huffman-revisited-part-4-multi-bytes.html)
#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![cfg_attr(docsrs, feature(doc_auto_cfg))]

mod compress;
mod decompress;
mod tables;

pub use compress::{compress_to_vec, Compressor, StoredOnlyCompressor};
pub use decompress::{
    decompress_to_vec, decompress_to_vec_bounded, BoundedDecompressionError, DecompressionError,
    Decompressor,
};

/// Build a length limited huffman tree.
///
/// Dynamic programming algorithm from fpnge.
#[doc(hidden)]
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

const fn compute_codes<const NSYMS: usize>(lengths: &[u8; NSYMS]) -> Option<[u16; NSYMS]> {
    let mut codes = [0u16; NSYMS];

    let mut code = 0u32;

    let mut len = 1;
    while len <= 16 {
        let mut i = 0;
        while i < lengths.len() {
            if lengths[i] == len {
                codes[i] = (code as u16).reverse_bits() >> (16 - len);
                code += 1;
            }
            i += 1;
        }
        code <<= 1;
        len += 1;
    }

    if code == 2 << 16 {
        Some(codes)
    } else {
        None
    }
}

// #[cfg(test)]
// mod tests {

//     #[test]
//     fn header() {
//         use crate::tables::HUFFMAN_LENGTHS;

//         let mut runs = Vec::new();
//         let mut run_start = None;
//         for i in 0..HUFFMAN_LENGTHS.len() {
//             if run_start.is_none() {
//                 run_start = Some(i);
//             } else if HUFFMAN_LENGTHS[i] != HUFFMAN_LENGTHS[i - 1] {
//                 runs.push((HUFFMAN_LENGTHS[i - 1], i - run_start.unwrap()));
//                 run_start = Some(i);
//             }
//         }
//         runs.push((
//             *HUFFMAN_LENGTHS.last().unwrap(),
//             HUFFMAN_LENGTHS.len() - run_start.unwrap(),
//         ));

//         let mut freqs = [0u64; 14];
//         freqs[1] = 1;

//         let mut symbols = Vec::new();

//         for (val, mut count) in runs {
//             if count < 4 {
//                 freqs[val as usize] += count as u64;
//                 for _ in 0..count {
//                     symbols.push((val, 0));
//                 }
//             } else {
//                 symbols.push((val, 0));
//                 freqs[val as usize] += 1;
//                 count -= 1;

//                 while count >= 9 {
//                     freqs[13] += 1;
//                     count -= 6;
//                     symbols.push((16, 3));
//                 }
//                 // match count {
//                 //     0 ..= 2 => unreachable!("Shouldn't happen"),
//                 //     3 ..= 6 => freqs[13] += 1,
//                 //     7 ..= 12 => freqs[13] += 2,
//                 //     13 .. => unreachable!("Shouldn't happen"),
//                 // }
//                 match count {
//                     0..=2 => unreachable!("Shouldn't happen"),
//                     3..=6 => symbols.push((16, count - 3)),
//                     7 => {
//                         symbols.push((16, 0));
//                         symbols.push((16, 1));
//                     }
//                     8 => {
//                         symbols.push((16, 1));
//                         symbols.push((16, 1));
//                     }
//                     9.. => unreachable!("Shouldn't happen"),
//                 }
//             }
//         }
//         println!("symbols = {:?}", symbols);

//         let mut lengths = [0u8; 13];
//         crate::compute_code_lengths(&freqs[1..], &[1; 13], &[7; 13], &mut lengths);

//         let mut bits_used = 0;
//         for i in 0..13 {
//             bits_used += u64::from(freqs[i + 1]) * u64::from(lengths[i]);
//         }
//         println!("Bits used: {}", bits_used);

//         println!("{:?}", lengths);

//         panic!();
//     }
// }
