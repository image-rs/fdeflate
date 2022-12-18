#![no_main]
use libfuzzer_sys::fuzz_target;
use miniz_oxide::inflate::TINFLStatus;
use std::io::{Cursor, Read};

fuzz_target!(|input: &[u8]| {
    if input.is_empty() {
        return;
    }

    match fdeflate::decompress_to_vec(&input) {
        Ok(decompressed) => {
            let decompressed2 = miniz_oxide::inflate::decompress_to_vec_zlib(&input).unwrap();
            assert_eq!(decompressed, decompressed2);
        }
        Err(fdeflate::DecompressionError::NotFDeflate) => {}
        Err(fdeflate::DecompressionError::InvalidHlit) => {
            // TODO: Remove once https://github.com/Frommi/miniz_oxide/issues/130 is fixed.
        }
        Err(err) => match miniz_oxide::inflate::decompress_to_vec_zlib(&input) {
            Err(r)
                if r.status == TINFLStatus::Failed
                    || r.status == TINFLStatus::FailedCannotMakeProgress => {}
            r => {
                panic!("fdeflate: {:?}, miniz_oxide: {:?}", err, r)
            }
        },
    }
});
