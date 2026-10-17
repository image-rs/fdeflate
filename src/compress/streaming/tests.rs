use super::*;

fn random(len: usize) -> Vec<u8> {
    let mut state = 0x12345678u32;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state as u8
        })
        .collect()
}

fn decode(data: &[u8], zlib: bool) -> Vec<u8> {
    if zlib {
        miniz_oxide::inflate::decompress_to_vec_zlib(data)
    } else {
        miniz_oxide::inflate::decompress_to_vec(data)
    }
    .unwrap()
}

fn compress(
    encoder: &mut StreamingCompressor,
    input: &[u8],
    input_chunk: usize,
    output_chunk: usize,
) -> Vec<u8> {
    let capacities = (encoder.input.data.capacity(), encoder.pending.capacity());
    let mut encoded = Vec::new();
    let mut output = vec![0; output_chunk];
    let mut position = 0;
    loop {
        let end = (position + input_chunk).min(input.len());
        let flush = if end == input.len() {
            CompressionFlush::Finish
        } else {
            CompressionFlush::None
        };
        let p = encoder.process(&input[position..end], &mut output, flush);
        position += p.consumed;
        encoded.extend_from_slice(&output[..p.produced]);
        assert_eq!(
            capacities,
            (encoder.input.data.capacity(), encoder.pending.capacity())
        );
        assert!(encoder.input.data.len() <= INPUT_CAPACITY);
        assert!(encoder.pending.len() <= OUTPUT_CAPACITY);
        if p.status.unwrap() == CompressionStatus::StreamEnd {
            break;
        }
        assert!(p.consumed != 0 || p.produced != 0);
    }
    assert_eq!(position, input.len());
    encoded
}

#[test]
fn cross_decoder_levels_chunking_and_history() {
    let mut input = random(BATCH_SIZE * 2 + 317);
    // Exercise distant references across automatic block boundaries.
    for i in 32768..input.len() {
        if i % 131 < 100 {
            input[i] = input[i - 32768];
        }
    }
    for zlib in [false, true] {
        for level in 0..=9 {
            let mut encoder = StreamingCompressor::new(level, zlib);
            for (ic, oc) in [(1, 4096), (4093, 1), (input.len(), input.len() + 1024)] {
                encoder.reset(zlib);
                let encoded = compress(&mut encoder, &input, ic, oc);
                assert_eq!(
                    decode(&encoded, zlib),
                    input,
                    "level={level} chunks={ic},{oc}"
                );
            }
        }
        let mut encoder = StreamingCompressor::new_rle(zlib);
        let input = vec![0; BATCH_SIZE * 3];
        assert_eq!(decode(&compress(&mut encoder, &input, 7, 1), zlib), input);
    }
}

#[test]
fn one_block_survives_window_slides_rebasing_and_empty_finish() {
    let mut input = random(1024 * 1024);
    for i in 8192..input.len() {
        input[i] = input[i - 8192];
    }
    for level in [1, 6, 9] {
        for zlib in [false, true] {
            let mut encoder = StreamingCompressor::new(level, zlib);
            let capacities = (encoder.input.data.capacity(), encoder.pending.capacity());
            let mut output = vec![0; OUTPUT_CAPACITY];
            let mut encoded = Vec::new();
            for chunk in input.chunks(BATCH_SIZE) {
                let p = encoder.process(chunk, &mut output, CompressionFlush::None);
                assert_eq!(p.consumed, chunk.len());
                assert_eq!(p.status.unwrap(), CompressionStatus::NeedsInput);
                encoded.extend_from_slice(&output[..p.produced]);
                // Tokens and saved literal bytes must survive rebasing too.
                encoder.rebase();
                assert_eq!(
                    capacities,
                    (encoder.input.data.capacity(), encoder.pending.capacity())
                );
            }
            assert_eq!(encoded.len(), if zlib { 2 } else { 0 });
            encoded.extend(compress(&mut encoder, &[], 1, 1));
            // The first DEFLATE block is final: sliding the input window did
            // not force any intervening Huffman blocks.
            assert_eq!(encoded[if zlib { 2 } else { 0 }] & 1, 1);
            assert_eq!(decode(&encoded, zlib), input);
        }
    }
}

#[test]
fn token_limit_applies_to_repetitive_input_without_a_flush() {
    let input = vec![b'A'; 8 * 1024 * 1024];
    for level in [1, 6, 9] {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut output = [0; 1];
        let p = encoder.process(&input, &mut output, CompressionFlush::None);
        assert_eq!(p.status.unwrap(), CompressionStatus::NeedsOutput);
        assert!(p.consumed > BATCH_SIZE && p.consumed < input.len());
        assert_eq!(p.produced, 1);
        assert_eq!(output[0] & 1, 0);
        let mut encoded = output.to_vec();
        encoded.extend(compress(&mut encoder, &input[p.consumed..], 65537, 127));
        assert_eq!(decode(&encoded, false), input);
    }
}

#[test]
fn literal_budget_splits_blocks_after_the_raw_prefix_has_left_the_window() {
    let mut input = random(1024 * 1024);
    let text = b"literal storage is bounded independently of match tokens\n";
    for (i, byte) in input.iter_mut().enumerate() {
        if i % 32768 >= 16384 {
            *byte = text[i % text.len()];
        }
    }
    for level in 1..=9 {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut output = [0; 1];
        let p = encoder.process(&input, &mut output, CompressionFlush::None);
        assert_eq!(p.status.unwrap(), CompressionStatus::NeedsOutput);
        assert!(p.consumed > BATCH_SIZE && p.consumed < input.len());
        assert_eq!(p.produced, 1);
        let mut encoded = output.to_vec();
        encoded.extend(compress(&mut encoder, &input[p.consumed..], 70003, 131));
        assert_eq!(decode(&encoded, false), input);
    }
}

#[test]
fn flushes_drain_saved_blocks_and_full_flush_restarts_history() {
    let input = b"saved literals and repeated matches across input batches\n".repeat(6000);
    for level in [1, 6, 9] {
        for flush in [
            CompressionFlush::Partial,
            CompressionFlush::Sync,
            CompressionFlush::Full,
        ] {
            let mut encoder = StreamingCompressor::new(level, false);
            let p = encoder.process(&input, &mut [0; 1], CompressionFlush::None);
            assert_eq!(p.consumed, input.len());
            assert_eq!(p.produced, 0);
            assert_eq!(p.status.unwrap(), CompressionStatus::NeedsInput);
            let mut encoded = Vec::new();
            loop {
                let mut output = [0; 1];
                let p = encoder.process(&[], &mut output, flush);
                encoded.extend_from_slice(&output[..p.produced]);
                if p.status.unwrap() == CompressionStatus::FlushComplete {
                    break;
                }
                assert_ne!(p.produced, 0);
            }
            let suffix = compress(&mut encoder, &input, 70003, 23);
            if flush == CompressionFlush::Full {
                assert_eq!(decode(&suffix, false), input);
            }
            encoded.extend(suffix);
            assert_eq!(decode(&encoded, false), input.repeat(2));
        }
    }
}

#[test]
fn incompressible_stream_is_bounded_and_emitted_before_finish() {
    let input = random(BATCH_SIZE);
    let mut encoder = StreamingCompressor::new(6, false);
    let capacities = (encoder.input.data.capacity(), encoder.pending.capacity());
    let mut output = vec![0; OUTPUT_CAPACITY];
    for _ in 0..128 {
        let p = encoder.process(&input, &mut output, CompressionFlush::None);
        assert_eq!(p.consumed, input.len());
        assert!(p.produced > 0);
        assert_eq!(p.status.unwrap(), CompressionStatus::NeedsInput);
        assert_eq!(
            capacities,
            (encoder.input.data.capacity(), encoder.pending.capacity())
        );
    }
}

#[test]
fn repeated_flushes_and_full_restart_with_one_byte_output() {
    let input = random(100000);
    for level in 0..=9 {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut encoded = Vec::new();
        let mut position = 0;
        loop {
            let mut output = [0; 1];
            let p = encoder.process(&input[position..], &mut output, CompressionFlush::Full);
            position += p.consumed;
            encoded.extend_from_slice(&output[..p.produced]);
            if p.status.unwrap() == CompressionStatus::FlushComplete {
                break;
            }
            assert!(p.consumed != 0 || p.produced != 0);
        }
        assert_eq!(position, input.len());
        for flush in [
            CompressionFlush::Partial,
            CompressionFlush::Sync,
            CompressionFlush::Full,
        ] {
            let p = encoder.process(&[], &mut [0; 1], flush);
            assert_eq!(p.status.unwrap(), CompressionStatus::FlushComplete);
            assert_eq!(p.produced, 0);
        }
        let suffix = compress(&mut encoder, &input, 101, 13);
        assert_eq!(decode(&suffix, false), input);
        encoded.extend_from_slice(&suffix);
        assert_eq!(decode(&encoded, false), input.repeat(2));
    }
}

#[test]
fn flush_request_survives_zero_output_and_none_retry() {
    for flush in [
        CompressionFlush::Partial,
        CompressionFlush::Sync,
        CompressionFlush::Full,
        CompressionFlush::Finish,
    ] {
        let mut encoder = StreamingCompressor::new(6, true);
        let mut output = vec![0; OUTPUT_CAPACITY];
        let first = encoder.process(b"hello", &mut output, CompressionFlush::None);
        let mut encoded = output[..first.produced].to_vec();
        assert_eq!(first.consumed, 5);
        assert_eq!(
            encoder.process(&[], &mut [], flush).status.unwrap(),
            CompressionStatus::NeedsOutput
        );
        loop {
            let p = encoder.process(&[], &mut output[..1], CompressionFlush::None);
            encoded.extend_from_slice(&output[..p.produced]);
            if p.status.unwrap() != CompressionStatus::NeedsOutput {
                break;
            }
        }
        if flush != CompressionFlush::Finish {
            encoded.extend(compress(&mut encoder, &[], 1, 1));
        }
        assert_eq!(decode(&encoded, true), b"hello");
    }
}

#[test]
fn frequent_rebasing_preserves_hash_chain_slots() {
    let mut input = random(BATCH_SIZE * 5);
    for i in 32768..input.len() {
        input[i] = input[i - 32768];
    }
    for level in 1..=9 {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut expected = StreamingCompressor::new(level, false);
        let mut encoded = Vec::new();
        let mut output = vec![0; OUTPUT_CAPACITY];
        for chunk in input.chunks(70003) {
            let p = encoder.process(chunk, &mut output, CompressionFlush::Sync);
            assert_eq!(p.consumed, chunk.len());
            assert_eq!(p.status.unwrap(), CompressionStatus::FlushComplete);
            let actual = output[..p.produced].to_vec();
            let q = expected.process(chunk, &mut output, CompressionFlush::Sync);
            assert_eq!(
                actual,
                output[..q.produced],
                "rebasing changed compression at level {level}"
            );
            encoded.extend(actual);
            encoder.rebase();
        }
        encoded.extend(compress(&mut encoder, &[], 1, 7));
        assert_eq!(decode(&encoded, false), input);
    }
}

#[test]
fn reset_discards_pending_output_and_retains_allocations() {
    let mut encoder = StreamingCompressor::new(6, true);
    let input = random(BATCH_SIZE * 2);
    encoder.process(&input, &mut [0; 3], CompressionFlush::Finish);
    let capacities = (encoder.input.data.capacity(), encoder.pending.capacity());
    for zlib in [false, true] {
        encoder.reset(zlib);
        assert_eq!(
            capacities,
            (encoder.input.data.capacity(), encoder.pending.capacity())
        );
        assert_eq!(
            decode(&compress(&mut encoder, b"new stream", 1, 1), zlib),
            b"new stream"
        );
        assert_eq!(
            encoder
                .process(b"extra", &mut [0; 1], CompressionFlush::None)
                .status,
            Err(CompressionError::Finished)
        );
    }
}

#[test]
#[ignore = "processes more than 4 GiB; run explicitly in release mode"]
fn long_stream_rebases_without_overflow() {
    let input = vec![b'A'; BATCH_SIZE];
    for level in [0, 1, 6] {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut output = vec![0; OUTPUT_CAPACITY];
        let mut decoder =
            miniz_oxide::inflate::stream::InflateState::new_boxed(miniz_oxide::DataFormat::Raw);
        let mut decoded = vec![0; BATCH_SIZE + 258];
        let mut total_decoded = 0u64;
        let mut verify = |bytes: &[u8], finished: bool| {
            let mut position = 0;
            loop {
                let result = miniz_oxide::inflate::stream::inflate(
                    &mut decoder,
                    &bytes[position..],
                    &mut decoded,
                    miniz_oxide::MZFlush::None,
                );
                position += result.bytes_consumed;
                assert!(decoded[..result.bytes_written].iter().all(|&b| b == b'A'));
                total_decoded += result.bytes_written as u64;
                if result.status == Ok(miniz_oxide::MZStatus::StreamEnd) {
                    assert!(finished);
                    break;
                }
                if result.bytes_consumed == 0 && result.bytes_written == 0 {
                    assert!(!finished);
                    assert_eq!(result.status, Err(miniz_oxide::MZError::Buf));
                    break;
                }
                if result.status == Err(miniz_oxide::MZError::Buf) {
                    // An empty-input drain can produce buffered output and
                    // report input starvation in the same miniz call.
                    assert!(!finished);
                    assert_eq!(position, bytes.len());
                } else {
                    assert_eq!(result.status, Ok(miniz_oxide::MZStatus::Ok));
                }
            }
            assert_eq!(position, bytes.len());
        };
        let mut count = 0u64;
        for _ in 0..(u32::MAX as u64 / BATCH_SIZE as u64 + 3) {
            let p = encoder.process(&input, &mut output, CompressionFlush::None);
            assert_eq!(p.consumed, input.len());
            assert_eq!(p.status.unwrap(), CompressionStatus::NeedsInput);
            count += p.consumed as u64;
            verify(&output[..p.produced], false);
        }
        assert!(count > u32::MAX as u64);
        let p = encoder.process(&[], &mut output, CompressionFlush::Finish);
        assert_eq!(p.status.unwrap(), CompressionStatus::StreamEnd);
        verify(&output[..p.produced], true);
        assert_eq!(total_decoded, count);
    }
}

#[test]
fn pending_flush_can_be_upgraded_without_losing_input() {
    let input = random(BATCH_SIZE);
    for level in [0, 1, 6, 9] {
        let mut encoder = StreamingCompressor::new(level, false);
        let mut encoded = Vec::new();
        let mut byte = [0; 1];
        let p = encoder.process(&input, &mut byte, CompressionFlush::Partial);
        assert_eq!(p.consumed, input.len());
        encoded.extend_from_slice(&byte[..p.produced]);
        encoder.process(&[], &mut [], CompressionFlush::Sync);
        encoder.process(&[], &mut [], CompressionFlush::Full);
        loop {
            let p = encoder.process(&[], &mut byte, CompressionFlush::None);
            encoded.extend_from_slice(&byte[..p.produced]);
            if p.status.unwrap() != CompressionStatus::NeedsOutput {
                break;
            }
        }
        let suffix = compress(&mut encoder, &input, 7, 1);
        assert_eq!(decode(&suffix, false), input);
        encoded.extend(suffix);
        assert_eq!(decode(&encoded, false), input.repeat(2));
    }
}

#[test]
fn every_flush_exposes_all_input_before_finishing() {
    for level in 0..=9 {
        for flush in [
            CompressionFlush::Partial,
            CompressionFlush::Sync,
            CompressionFlush::Full,
        ] {
            for size in [1, 258, BATCH_SIZE - 1, BATCH_SIZE + 1] {
                let input = random(size);
                let mut encoder = StreamingCompressor::new(level, false);
                let mut encoded = Vec::new();
                let mut position = 0;
                loop {
                    let mut byte = [0; 1];
                    let p = encoder.process(&input[position..], &mut byte, flush);
                    position += p.consumed;
                    encoded.extend_from_slice(&byte[..p.produced]);
                    if p.status.unwrap() == CompressionStatus::FlushComplete {
                        break;
                    }
                    assert!(p.consumed != 0 || p.produced != 0);
                }
                let mut decoder = miniz_oxide::inflate::stream::InflateState::new_boxed(
                    miniz_oxide::DataFormat::Raw,
                );
                let mut decoded = vec![0; input.len() + 1];
                let (mut consumed, mut produced) = (0, 0);
                loop {
                    let result = miniz_oxide::inflate::stream::inflate(
                        &mut decoder,
                        &encoded[consumed..],
                        &mut decoded[produced..],
                        miniz_oxide::MZFlush::None,
                    );
                    consumed += result.bytes_consumed;
                    produced += result.bytes_written;
                    if result.bytes_consumed == 0 && result.bytes_written == 0 {
                        assert_eq!(result.status, Err(miniz_oxide::MZError::Buf));
                        break;
                    }
                    assert_eq!(result.status, Ok(miniz_oxide::MZStatus::Ok));
                }
                assert_eq!(position, input.len());
                assert_eq!(produced, input.len(), "level={level} flush={flush:?}");
                assert_eq!(&decoded[..produced], input);
            }
        }
    }
}
