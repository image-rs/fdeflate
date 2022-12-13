#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::{Cursor, Read};

fuzz_target!(|input: &[u8]| {
    let mut decompressed = vec![0; input.len() * 2 + 2];
    let mut decompressor = fdeflate::Decompressor::new();

    match decompressor.read(&input, &mut decompressed, true) {
        Ok(n) => {
            assert!(decompressor.done());
            let mut decompressed2 = Vec::new();
            flate2::bufread::ZlibDecoder::new(Cursor::new(&input))
                .read_to_end(&mut decompressed2)
                .unwrap();
            assert_eq!(decompressed[..n.1], decompressed2);
        }
        Err(fdeflate::DecompressionError::CorruptData) => {
            match flate2::Decompress::new(true).decompress_vec(
                &input,
                &mut Vec::with_capacity(input.len() * 2 + 2),
                flate2::FlushDecompress::Finish,
            ) {
                Ok(flate2::Status::BufError) | Err(_) => {}
                r => panic!("{:?}", r),
            }
        }
        Err(fdeflate::DecompressionError::NotFDeflate) => {}
    }
});
