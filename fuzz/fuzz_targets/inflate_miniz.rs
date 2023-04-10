#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate miniz_oxide;

fuzz_target!(|input: (u8, Vec<u8>)| {
    let compression_level = input.0;
    let data = input.1;
    let compressed = miniz_oxide::deflate::compress_to_vec_zlib(&data, compression_level);
    let decompressed = fdeflate::decompress_to_vec(&compressed).expect("Decompression failed!");
    assert_eq!(data, decompressed);
});
