//! This fuzz target tests that reading bytes out of the decompressor one at a time always produces
//! valid output.

#![no_main]
#[macro_use]
extern crate libfuzzer_sys;
extern crate miniz_oxide;

fuzz_target!(|input: (u8, Vec<u8>)| {
    let compression_level = input.0;
    let data = input.1;
    if data.is_empty() {
        return;
    }
    let compressed = miniz_oxide::deflate::compress_to_vec_zlib(&data, compression_level);

    let mut decompressed = Vec::new();
    let mut input_index = 0;
    let mut decoder = fdeflate::Decompressor::new();
    loop {
        let output_position = decompressed.len();
        if output_position < data.len() {
            decompressed.push(1);
        }

        let (consumed, produced) = decoder
            .read(
                &compressed[input_index..],
                &mut decompressed,
                output_position,
            )
            .expect("Decompression failed!");

        input_index += consumed;
        assert_eq!(produced, 1);
        if decoder.is_done() {
            break;
        }
    }

    assert_eq!(data, decompressed);
});
