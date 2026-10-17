#![forbid(unsafe_code)]

use fdeflate::{DecompressionError, DecompressionStatus, Decompressor, Format};

#[test]
fn checksum_error_reports_progress_and_requires_reset() {
    let input = vec![42; 20000];
    let mut encoded = miniz_oxide::deflate::compress_to_vec_zlib(&input, 6);
    *encoded.last_mut().unwrap() ^= 1;
    let mut decoder = Decompressor::new();
    let mut output = vec![0; input.len()];
    let result = decoder.read_progress(&encoded, &mut output, 0);
    assert_eq!(result.consumed, encoded.len());
    assert_eq!(result.produced, input.len());
    assert!(matches!(
        result.status,
        Err(DecompressionError::WrongChecksum)
    ));
    assert_eq!(output, input);
    let result = decoder.read_progress(&encoded, &mut output, 0);
    assert_eq!((result.consumed, result.produced), (0, 0));
    assert!(matches!(
        result.status,
        Err(DecompressionError::WrongChecksum)
    ));
    decoder.reset(Format::Zlib);
    *encoded.last_mut().unwrap() ^= 1;
    let result = decoder.read_progress(&encoded, &mut output, 0);
    assert_eq!(result.status.unwrap(), DecompressionStatus::Done);
    assert_eq!(result.produced, input.len());
}

#[test]
fn error_after_valid_block_reports_its_output() {
    // Non-final stored block with one byte, then an invalid block type.
    let encoded = [0, 1, 0, 0xfe, 0xff, b'x', 7, 0];
    let mut decoder = Decompressor::new_with_format(Format::Raw);
    let mut output = [0; 100];
    let result = decoder.read_progress(&encoded, &mut output, 0);
    assert!(matches!(
        result.status,
        Err(DecompressionError::InvalidBlockType)
    ));
    assert_eq!(result.consumed, encoded.len());
    assert_eq!(result.produced, 1);
    assert_eq!(output[0], b'x');
}

#[test]
fn empty_output_can_consume_empty_blocks_and_checksum() {
    for format in [Format::Raw, Format::Zlib] {
        let encoded = if format == Format::Raw {
            miniz_oxide::deflate::compress_to_vec(&[], 0)
        } else {
            miniz_oxide::deflate::compress_to_vec_zlib(&[], 0)
        };
        let mut decoder = Decompressor::new_with_format(format);
        for (index, byte) in encoded.iter().enumerate() {
            let result = decoder.read_progress(&[*byte], &mut [], 0);
            assert_eq!((result.consumed, result.produced), (1, 0));
            assert_eq!(
                result.status.unwrap(),
                if index + 1 == encoded.len() {
                    DecompressionStatus::Done
                } else {
                    DecompressionStatus::NeedsInput
                }
            );
        }
    }
}

#[test]
fn reserved_fixed_symbol_is_error_with_committed_progress() {
    let input = [75, 28, 5, 163, 96, 20, 140, 221, 0, 0];
    let mut decoder = Decompressor::new_with_format(Format::Raw);
    let mut output = [0; 2048];
    let result = decoder.read_progress(&input, &mut output, 0);
    assert!(matches!(
        result.status,
        Err(DecompressionError::InvalidLiteralLengthCode)
    ));
    assert!(result.produced > 0);
    assert!(output[..result.produced].iter().all(|&x| x == b'a'));
}

#[test]
fn both_reserved_fixed_symbols_fail_in_fast_and_slow_paths() {
    for symbol in [286u16, 287] {
        let literal = (0x30u16 + u16::from(b'a')).reverse_bits() >> 8;
        let reserved = (0xc0u16 + symbol - 280).reverse_bits() >> 8;
        let bits = 3u32 | (u32::from(literal) << 3) | (u32::from(reserved) << 11);
        for padding in [0, 16] {
            let mut input = bits.to_le_bytes()[..3].to_vec();
            input.resize(input.len() + padding, 0);
            for chunk in [1, input.len()] {
                let mut decoder = Decompressor::new_with_format(Format::Raw);
                let mut output = [0; 256];
                let mut produced = 0;
                let mut errored = false;
                for part in input.chunks(chunk) {
                    let result = decoder.read_progress(part, &mut output, produced);
                    produced += result.produced;
                    if let Err(error) = result.status {
                        assert_eq!(error, DecompressionError::InvalidLiteralLengthCode);
                        errored = true;
                        break;
                    }
                }
                assert!(errored);
                assert_eq!(produced, 1);
                assert_eq!(output[0], b'a');
            }
        }
    }
}
