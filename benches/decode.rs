#![feature(test)]

extern crate test;

use fdeflate::{Decompressor, Format};
use test::{black_box, Bencher};

fn data(kind: &str, len: usize) -> Vec<u8> {
    let mut state = 0x12345678u32;
    (0..len)
        .map(|i| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            match kind {
                "literals" => (state % 16) as u8,
                "backrefs" => b"abcde"[i % 5],
                "stored" => state as u8,
                _ => {
                    if i % 512 < 256 {
                        (state % 32) as u8
                    } else {
                        b"abcde"[i % 5]
                    }
                }
            }
        })
        .collect()
}

fn decode(
    b: &mut Bencher,
    kind: &str,
    format: Format,
    input_chunk: usize,
    output_chunk: usize,
    len: usize,
) {
    let data = data(kind, len);
    let level = if kind == "stored" { 0 } else { 6 };
    let compressed = if format == Format::Raw {
        miniz_oxide::deflate::compress_to_vec(&data, level)
    } else {
        miniz_oxide::deflate::compress_to_vec_zlib(&data, level)
    };
    let mut decoder = Decompressor::new_with_format(format);
    // Spare space also lets the old implementation finish, for before/after comparisons.
    let mut output = vec![0; data.len() + 64];
    b.bytes = len as u64;
    b.iter(|| {
        decoder.reset(format);
        let (mut ip, mut op) = (0usize, 0usize);
        while !decoder.is_done() {
            let input_end = ip.saturating_add(input_chunk).min(compressed.len());
            let output_end = op.saturating_add(output_chunk).min(output.len());
            let (consumed, produced) = decoder
                .read(
                    black_box(&compressed[ip..input_end]),
                    &mut output[..output_end],
                    op,
                )
                .unwrap();
            assert!(consumed != 0 || produced != 0 || decoder.is_done());
            ip += consumed;
            op += produced;
        }
        assert_eq!(ip, compressed.len());
        assert_eq!(op, data.len());
        black_box(&output[..op]);
    });
    assert_eq!(&output[..data.len()], data);
}

macro_rules! bench {
    ($name:ident, $kind:literal, $format:ident, $input:expr, $output:expr, $len:expr) => {
        #[bench]
        fn $name(b: &mut Bencher) {
            decode(b, $kind, Format::$format, $input, $output, $len);
        }
    };
}

bench!(
    raw_literals,
    "literals",
    Raw,
    usize::MAX,
    usize::MAX,
    1 << 20
);
bench!(
    zlib_literals,
    "literals",
    Zlib,
    usize::MAX,
    usize::MAX,
    1 << 20
);
bench!(
    raw_backrefs,
    "backrefs",
    Raw,
    usize::MAX,
    usize::MAX,
    1 << 20
);
bench!(
    zlib_backrefs,
    "backrefs",
    Zlib,
    usize::MAX,
    usize::MAX,
    1 << 20
);
bench!(raw_stored, "stored", Raw, usize::MAX, usize::MAX, 1 << 20);
bench!(zlib_stored, "stored", Zlib, usize::MAX, usize::MAX, 1 << 20);
bench!(raw_mixed, "mixed", Raw, usize::MAX, usize::MAX, 1 << 20);
bench!(zlib_mixed, "mixed", Zlib, usize::MAX, usize::MAX, 1 << 20);
bench!(raw_chunks, "mixed", Raw, 65536, 65536, 1 << 20);
bench!(zlib_chunks, "mixed", Zlib, 65536, 65536, 1 << 20);
bench!(bytewise_input, "mixed", Zlib, 1, usize::MAX, 16384);
bench!(bytewise_output, "mixed", Zlib, usize::MAX, 1, 16384);
bench!(small_reset_raw, "backrefs", Raw, usize::MAX, usize::MAX, 32);
bench!(
    small_reset_zlib,
    "backrefs",
    Zlib,
    usize::MAX,
    usize::MAX,
    32
);
