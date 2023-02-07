#![no_main]
use libfuzzer_sys::fuzz_target;
use miniz_oxide::inflate::TINFLStatus;
use std::io::Read;

fuzz_target!(|input: &[u8]| {
    if input.is_empty() {
        return;
    }

    match fdeflate::decompress_to_vec(&input) {
        Ok(decompressed) => {
            if decompressed.is_empty() {
                return;
            }
            let Ok(decompressed2) = miniz_oxide::inflate::decompress_to_vec_zlib(&input) else {return};
            // let min_len = decompressed.len().min(decompressed2.len());
            // assert_eq!(decompressed[..min_len], decompressed2[..min_len]);
            // assert_eq!(decompressed.len(), decompressed2.len());
            assert_eq!(decompressed, decompressed2);
        }
        Err(fdeflate::DecompressionError::BadLiteralLengthHuffmanTree) => {}
        Err(err) => match miniz_oxide::inflate::decompress_to_vec_zlib(&input) {
            Err(r)
                if r.status == TINFLStatus::Failed
                    || r.status == TINFLStatus::FailedCannotMakeProgress => {}
            r => {
                if let Err(e) = flate2::read::ZlibDecoder::new(std::io::Cursor::new(input))
                    .read_to_end(&mut Vec::new())
                {
                    panic!("fdeflate: {:?}, miniz_oxide: {:?}, zlib: {:?}", err, r, e);
                }
            }
        },
    }
});
