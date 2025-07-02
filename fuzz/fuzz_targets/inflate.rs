#![no_main]
use libfuzzer_sys::fuzz_target;
use miniz_oxide::inflate::TINFLStatus;

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
            assert_eq!(decompressed, decompressed2);
        }
        Err(fdeflate::DecompressionError::BadCodeLengthHuffmanTree) => {}
        Err(fdeflate::DecompressionError::BadLiteralLengthHuffmanTree) => {}
        Err(fdeflate::DecompressionError::BadDistanceHuffmanTree) => {}
        Err(fdeflate::DecompressionError::InvalidDistanceCode) => {}
        Err(err) => match miniz_oxide::inflate::decompress_to_vec_zlib(&input) {
            Err(r)
                if r.status == TINFLStatus::Failed
                    || r.status == TINFLStatus::FailedCannotMakeProgress => {}
            r => {
                panic!("fdeflate: {:?}, miniz_oxide: {:?}", err, r);
            }
        },
    }
});
