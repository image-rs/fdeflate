#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::{Cursor, Read};

fuzz_target!(|data: Vec<Vec<u8>>| {
    let mut compressor = fdeflate::StoredOnlyCompressor::new(Cursor::new(Vec::new())).unwrap();
    for chunk in &data {
        compressor.write_data(&*chunk).unwrap();
    }
    let compressed = compressor.finish().unwrap().into_inner();

    let mut decompressed = Vec::new();
    flate2::bufread::ZlibDecoder::new(Cursor::new(&compressed))
        .read_to_end(&mut decompressed)
        .unwrap();
    assert_eq!(decompressed, data.concat());

    assert_eq!(
        fdeflate::StoredOnlyCompressor::<()>::compressed_size(decompressed.len()),
        compressed.len(),
    );
});
