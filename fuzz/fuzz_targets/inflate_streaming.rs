//! Exercise input ownership and forward progress with independent input/output schedules.
#![no_main]

use libfuzzer_sys::fuzz_target;

#[path = "../../src/decompress/tests/streaming_utils.rs"]
mod streaming_utils;
use streaming_utils::stream;

fuzz_target!(|input: (u8, Vec<u8>, Vec<u8>, Vec<u8>)| {
    let (control, mut data, input_schedule, output_schedule) = input;
    data.truncate(8192);
    let schedule = |bytes: &[u8]| {
        let mut chunks: Vec<usize> = bytes.iter().take(16).map(|&n| n as usize).collect();
        // Guarantee progress even for an empty or all-zero schedule.
        chunks.push(1);
        chunks
    };
    let input_chunks = schedule(&input_schedule);
    let output_chunks = schedule(&output_schedule);
    let format = if control & 1 == 0 {
        fdeflate::Format::Raw
    } else {
        fdeflate::Format::Zlib
    };
    let compressed = if format == fdeflate::Format::Raw {
        miniz_oxide::deflate::compress_to_vec(&data, (control >> 1) % 10)
    } else {
        miniz_oxide::deflate::compress_to_vec_zlib(&data, (control >> 1) % 10)
    };
    let reference = if format == fdeflate::Format::Raw {
        miniz_oxide::inflate::decompress_to_vec(&compressed)
    } else {
        miniz_oxide::inflate::decompress_to_vec_zlib(&compressed)
    }
    .unwrap();
    assert_eq!(reference, data);

    let mut with_trailing = compressed.clone();
    with_trailing.extend_from_slice(&compressed);
    with_trailing.extend_from_slice(&input_schedule);
    for chunks in [&[usize::MAX][..], &input_chunks] {
        let result = stream(
            &with_trailing,
            format,
            data.len(),
            chunks,
            &output_chunks,
            Some(compressed.len()),
        )
        .unwrap();
        assert!(result.done);
        assert_eq!(result.consumed, compressed.len());
        assert_eq!(result.output, reference);
    }
    let cut = control as usize % compressed.len();
    if let Ok(result) = stream(
        &compressed[..cut],
        format,
        data.len(),
        &input_chunks,
        &output_chunks,
        None,
    ) {
        assert!(!result.done, "accepted truncated stream");
    }

    // Arbitrary malformed input also exercises the parser, with bounded memory and calls.
    // Error variants can differ with lookahead, but successful completion must agree.
    let whole = stream(&data, format, 65536, &[usize::MAX], &[usize::MAX], None);
    let chunked = stream(&data, format, 65536, &input_chunks, &output_chunks, None);
    match (whole, chunked) {
        (Ok(a), Ok(b)) => {
            assert_eq!(a.done, b.done);
            if a.done {
                assert_eq!(a.output, b.output);
                assert_eq!(a.consumed, b.consumed);
            }
        }
        (Ok(a), Err(_)) => assert!(!a.done),
        (Err(_), Ok(b)) => assert!(!b.done),
        (Err(_), Err(_)) => (),
    }
});
