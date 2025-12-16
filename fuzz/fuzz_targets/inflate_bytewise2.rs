//! This fuzz target tests that feeding bytes into the decompressor one at a time always produces
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
    let mut output_position = 0;
    let mut decoder = fdeflate::Decompressor::new();
    loop {
        decompressed.resize(output_position + 1024, 0);

        let (consumed, produced) = decoder
            .read(
                &compressed[input_index..input_index + 1],
                &mut decompressed,
                output_position,
            )
            .expect("Decompression failed!");

        input_index += consumed;
        output_position += produced;
        if decoder.is_done() {
            break;
        }
    }

    decompressed.resize(output_position, 0);
    assert_eq!(data, decompressed);
});
