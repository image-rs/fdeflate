#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::{Cursor, Read};

fuzz_target!(|input: (&[u8], bool)| {
    let mut decompressed = vec![0; input.0.len() * 2 + 2];
    let mut decompressor = fdeflate::Decompressor::new();

    match decompressor.read(&input.0, &mut decompressed, input.1) {
        Ok(n) => {
            if decompressor.done() {
                let mut decompressed2 = Vec::new();
                flate2::bufread::ZlibDecoder::new(Cursor::new(&input.0))
                    .read_to_end(&mut decompressed2)
                    .unwrap();
                assert_eq!(decompressed[..n.1], decompressed2);
            } else {
                assert!(!input.1);
            }
        }
        Err(fdeflate::DecompressionError::CorruptData) => {
            match flate2::Decompress::new(true).decompress_vec(
                &input.0,
                &mut Vec::with_capacity(input.0.len() * 2 + 2),
                flate2::FlushDecompress::Finish,
            ) {
                Ok(flate2::Status::BufError) | Err(_) => {}
                r => panic!("{:?}", r),
            }
        }
        Err(fdeflate::DecompressionError::NotFDeflate) => {}
    }
});
