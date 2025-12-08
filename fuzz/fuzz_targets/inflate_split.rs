//! This fuzz target tests that decompressing a byte slice always produces the same output as
//! decompressing in two steps.

#![no_main]
#[macro_use]
extern crate libfuzzer_sys;
extern crate miniz_oxide;

fuzz_target!(|input: (Vec<u8>, Vec<u8>)| {
    let joined = input.0.iter().chain(&input.1).copied().collect::<Vec<u8>>();
    let full_output = fdeflate::decompress_to_vec(&joined);
    // println!("-------------");

    let split_output: Result<_, fdeflate::DecompressionError> = (|| {
        let mut decoder = fdeflate::Decompressor::new();
        let mut output = vec![0; 1024];
        let mut input_index = 0;
        let mut output_index = 0;
        while !decoder.is_done() && input_index < input.0.len() {
            // println!("input_index: {} {}", input_index, input.0.len());
            let (consumed, produced) =
                decoder.read(&input.0[input_index..], &mut output, output_index)?;
            input_index += consumed;
            output_index += produced;
            if output_index == output.len() {
                output.resize(output_index + 32 * 1024, 0);
            }
            assert!(consumed > 0 || produced > 0 || decoder.is_done());
        }
        while !decoder.is_done() {
            let (consumed, produced) = decoder.read(
                &input.1[input_index - input.0.len()..],
                &mut output,
                output_index,
            )?;

            if !decoder.is_done() && consumed == 0 && produced == 0 {
                return Err(fdeflate::DecompressionError::InsufficientInput);
            }

            input_index += consumed;
            output_index += produced;
            output.resize(output_index + 32 * 1024, 0);
        }
        output.resize(output_index, 0);
        Ok(output)
    })();

    assert_eq!(full_output, split_output);
    //assert!(full_output == split_output || split_output.is_err() && full_output.is_err());
});
