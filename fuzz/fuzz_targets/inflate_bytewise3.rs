//! This fuzz target tests that feeding bytes into the decompressor one at a time always produces
//! valid output.

#![no_main]
#[macro_use]
extern crate libfuzzer_sys;
extern crate miniz_oxide;

#[path = "../../src/decompress/tests/test_utils.rs"]
mod test_utils;
use test_utils::decompress_by_chunks;

fuzz_target!(|input: &[u8]| {
    let r_whole = decompress_by_chunks(input, std::iter::repeat(input.len()));
    let r_bytewise = decompress_by_chunks(input, std::iter::repeat(1));
    match (r_whole, r_bytewise) {
        (Ok(output_whole), Ok(output_bytewise)) => assert_eq!(output_whole, output_bytewise),
        (Err(_e1), Err(_e2)) => (),
        (Ok(_), Err(e)) => panic!("Only byte-by-byte returned an error: {:?}", e),
        (Err(e), Ok(_)) => panic!("Only consume-whole returned an error: {:?}", e),
    }
});
