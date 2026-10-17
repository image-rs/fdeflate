use super::super::*;

#[path = "streaming_utils.rs"]
mod streaming_utils;
use streaming_utils::stream;

#[derive(Default)]
struct Bits {
    bytes: Vec<u8>,
    len: usize,
}

impl Bits {
    fn put(&mut self, value: u16, count: usize) {
        for i in 0..count {
            if self.len.is_multiple_of(8) {
                self.bytes.push(0);
            }
            self.bytes[self.len / 8] |= (((value >> i) & 1) as u8) << (self.len % 8);
            self.len += 1;
        }
    }
}

fn wrap(raw: &[u8], data: &[u8], format: Format) -> Vec<u8> {
    if format == Format::Raw {
        return raw.to_vec();
    }
    let mut input = vec![0x78, 0x01];
    input.extend_from_slice(raw);
    input.extend_from_slice(&simd_adler32::adler32(&data).to_be_bytes());
    input
}

// A complete comb-shaped literal tree: lengths 1, 2, ..., eof_bits - 1, followed
// by two codes of length eof_bits. Literal zero has code zero; EOF is all ones.
fn dynamic(eof_bits: usize, zeros: usize) -> Vec<u8> {
    dynamic_literals(eof_bits, &vec![0; zeros])
}

fn dynamic_literals(eof_bits: usize, literals: &[u8]) -> Vec<u8> {
    let mut lengths = [0u16; 258];
    for (i, length) in lengths.iter_mut().take(eof_bits - 1).enumerate() {
        *length = (i + 1) as u16;
    }
    lengths[eof_bits - 1] = eof_bits as u16;
    lengths[256] = eof_bits as u16;
    let mut bits = Bits::default();
    bits.put(0b101, 3); // Final dynamic block.
    bits.put(0, 5); // 257 literal/length codes.
    bits.put(0, 5); // One unused distance code.
    bits.put(15, 4); // All 19 code-length code lengths follow.
    for symbol in CLCL_ORDER {
        bits.put(if symbol < 16 { 4 } else { 0 }, 3);
    }
    for length in lengths {
        bits.put(length.reverse_bits() >> 12, 4);
    }
    for &literal in literals {
        let literal = literal as usize;
        assert!(literal < eof_bits);
        bits.put((1 << literal) - 1, (literal + 1).min(eof_bits));
    }
    bits.put((1 << eof_bits) - 1, eof_bits);
    bits.bytes
}

fn empty_fixed_blocks(count: usize) -> Vec<u8> {
    let mut bits = Bits::default();
    for _ in 0..count {
        bits.put(0b010, 3);
        bits.put(0, 7);
    }
    bits.put(0b011, 3);
    bits.put((0x30u16 + b'A' as u16).reverse_bits() >> 8, 8);
    bits.put(0, 7);
    bits.bytes
}

fn fixtures() -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut stored = vec![0, 0, 0, 0xff, 0xff]; // Empty nonfinal stored block.
    stored.extend_from_slice(&[1, 3, 0, 0xfc, 0xff, b'a', b'b', b'c']);
    let mut trailing_fixed = Bits::default();
    trailing_fixed.put(0b010, 3);
    trailing_fixed.put((0x30u16 + b'A' as u16).reverse_bits() >> 8, 8);
    trailing_fixed.put(0, 7);
    for i in 0..9 {
        trailing_fixed.put(if i == 8 { 0b011 } else { 0b010 }, 3);
        trailing_fixed.put(0, 7);
    }
    let mut trailing_stored = vec![0, 1, 0, 0xfe, 0xff, b'A'];
    trailing_stored.extend_from_slice(&[0, 0, 0, 0xff, 0xff, 1, 0, 0, 0xff, 0xff]);
    let mut trailing_dynamic = vec![0, 1, 0, 0xfe, 0xff, b'A'];
    trailing_dynamic.extend_from_slice(&dynamic(15, 0));
    vec![
        (vec![3, 0], vec![]),
        (vec![0x63, 0, 0], vec![0]),
        (empty_fixed_blocks(9), vec![b'A']),
        (stored, b"abc".to_vec()),
        (trailing_fixed.bytes, vec![b'A']),
        (trailing_stored, vec![b'A']),
        (trailing_dynamic, vec![b'A']),
        (dynamic(1, 9), vec![0; 9]),
        (dynamic(15, 9), vec![0; 9]),
        (dynamic_literals(15, &[14, 14]), vec![14, 14]),
    ]
}

#[test]
fn forward_progress_with_exact_output() {
    for input in [b"teststring".as_slice(), b"", b"\0"] {
        let zlib = crate::compress_to_vec(input);
        for format in [Format::Zlib, Format::Raw] {
            let compressed = if format == Format::Zlib {
                &zlib[..]
            } else {
                &zlib[2..zlib.len() - 4]
            };
            let mut decoder = Decompressor::new_with_format(format);
            let mut output = vec![0; input.len()];
            let (mut input_index, mut output_index) = (0, 0);
            while !decoder.is_done() {
                assert!(
                    input_index < compressed.len(),
                    "failed to finish with all input"
                );
                let (consumed, produced) = decoder
                    .read(&compressed[input_index..][..1], &mut output, output_index)
                    .unwrap();
                assert!(
                    consumed > 0,
                    "bytewise input stalled with exact-sized output"
                );
                input_index += consumed;
                output_index += produced;
            }
            assert_eq!(input_index, compressed.len());
            assert_eq!(output_index, input.len());
            assert_eq!(output, input);
        }
    }
}

#[test]
fn every_input_split_with_output_backpressure() {
    for (raw, data) in fixtures() {
        assert_eq!(miniz_oxide::inflate::decompress_to_vec(&raw).unwrap(), data);
        for format in [Format::Raw, Format::Zlib] {
            let compressed = wrap(&raw, &data, format);
            let end = compressed.len();
            let mut with_trailing = compressed.clone();
            with_trailing.extend_from_slice(&compressed); // Next stream must remain intact.
            with_trailing.extend_from_slice(b"trailing sentinel");
            for split in 0..=end {
                for output_chunks in [&[1][..], &[0, 1, 0, 3], &[usize::MAX]] {
                    for spare in [0, 1] {
                        let result = stream(
                            &with_trailing,
                            format,
                            data.len() + spare,
                            &[split, usize::MAX],
                            output_chunks,
                            Some(end),
                        )
                        .unwrap();
                        assert!(
                            result.done,
                            "format={format:?}, split={split}, spare={spare}"
                        );
                        assert_eq!(result.output, data);
                        assert_eq!(&with_trailing[result.consumed..], &with_trailing[end..]);
                    }
                }
            }
        }
    }
}

#[test]
fn all_end_code_lengths_and_alignments() {
    for eof_bits in 1..=15 {
        for zeros in 0..8 {
            let raw = dynamic(eof_bits, zeros);
            let data = vec![0; zeros];
            assert_eq!(miniz_oxide::inflate::decompress_to_vec(&raw).unwrap(), data);
            for format in [Format::Raw, Format::Zlib] {
                let compressed = wrap(&raw, &data, format);
                for input_chunks in [&[1][..], &[usize::MAX]] {
                    let result = stream(
                        &compressed,
                        format,
                        zeros,
                        input_chunks,
                        &[1],
                        Some(compressed.len()),
                    )
                    .unwrap();
                    assert!(
                        result.done,
                        "EOF bits={eof_bits}, zeros={zeros}, format={format:?}"
                    );
                    assert_eq!(result.output, data);
                }
            }
        }
    }
}

#[test]
fn partial_end_prefix_can_be_a_literal() {
    // This literal shares fourteen bits with EOF. With no output space the decoder must
    // retain those bits across input-starved calls, then return only the new byte when its
    // fifteenth bit proves that output is needed. Rewinding the old bytes would lose input.
    let raw = dynamic_literals(15, &[14]);
    for format in [Format::Raw, Format::Zlib] {
        let mut input = wrap(&raw, &[14], format);
        let end = input.len();
        input.extend_from_slice(b"trailing sentinel");
        let mut decoder = Decompressor::new_with_format(format);
        let mut pos = 0;
        loop {
            assert!(pos < end);
            let (consumed, produced) = decoder.read(&input[pos..][..1], &mut [], 0).unwrap();
            assert_eq!(produced, 0);
            pos += consumed;
            if consumed == 0 {
                break;
            }
        }
        assert_eq!(decoder.bits.nbits, 14);
        let mut output = [0];
        assert_eq!(
            decoder.read(&input[pos..], &mut output, 0),
            Ok((end - pos, 1))
        );
        assert!(decoder.is_done());
        assert_eq!(output, [14]);
    }
}

#[test]
fn queued_output_and_stored_copies_preserve_trailers() {
    for data in [
        vec![42; 2048],
        b"abcde".repeat(1000),
        (0..65536).map(|i| i as u8).collect(),
    ] {
        for level in [0, 6] {
            let raw = miniz_oxide::deflate::compress_to_vec(&data, level);
            for format in [Format::Raw, Format::Zlib] {
                let mut compressed = wrap(&raw, &data, format);
                let end = compressed.len();
                compressed.extend_from_slice(b"trailing sentinel");
                for (input_chunks, output_chunks) in [
                    (&[1][..], &[1][..]),
                    (&[3, 7, 64], &[0, 1, 17]),
                    (&[usize::MAX], &[1]),
                ] {
                    let result = stream(
                        &compressed,
                        format,
                        data.len(),
                        input_chunks,
                        output_chunks,
                        Some(end),
                    )
                    .unwrap();
                    assert!(result.done);
                    assert_eq!(result.output, data);
                }
            }
        }
    }
}

#[test]
fn truncated_and_invalid_streams_are_bounded() {
    for (raw, data) in fixtures() {
        for format in [Format::Raw, Format::Zlib] {
            let compressed = wrap(&raw, &data, format);
            for end in 0..compressed.len() {
                for chunks in [&[1][..], &[usize::MAX]] {
                    let result = stream(&compressed[..end], format, data.len(), chunks, &[1], None);
                    assert!(!result.is_ok_and(|r| r.done), "accepted truncated stream");
                }
            }
        }
    }
    for raw in [&[7, 0][..], &[1, 1, 0, 0, 0, 42]] {
        for chunks in [&[1][..], &[usize::MAX]] {
            assert!(stream(raw, Format::Raw, 1024, chunks, &[1], None).is_err());
        }
    }
}

#[test]
fn checksum_after_queued_output_and_format_resets() {
    let data = b"abcde".repeat(100);
    let raw = miniz_oxide::deflate::compress_to_vec(&data, 6);
    let mut zlib = wrap(&raw, &data, Format::Zlib);
    *zlib.last_mut().unwrap() ^= 1;
    for chunks in [&[1][..], &[usize::MAX]] {
        assert_eq!(
            stream(&zlib, Format::Zlib, data.len(), chunks, &[1], None).unwrap_err(),
            DecompressionError::WrongChecksum
        );
    }
    let mut decoder = Decompressor::new();
    for format in [Format::Raw, Format::Zlib, Format::Raw, Format::Zlib] {
        decoder.reset(format);
        let input = wrap(&raw, &data, format);
        let mut output = vec![0; data.len()];
        assert_eq!(
            decoder.read(&input, &mut output, 0).unwrap(),
            (input.len(), data.len())
        );
        assert!(decoder.is_done());
        assert_eq!(output, data);
        decoder.ignore_adler32(); // Reset must restore validation for the next zlib stream.
    }
    decoder.reset(Format::Zlib);
    assert_eq!(
        decoder.read(&zlib, &mut vec![0; data.len()], 0),
        Err(DecompressionError::WrongChecksum)
    );
}

#[test]
fn reset_discards_partial_headers_and_queued_output() {
    let compressed = crate::compress_to_vec_rle(&[42; 2048]);
    let mut decoder = Decompressor::new();
    for format in [Format::Raw, Format::Zlib] {
        decoder.reset(Format::Zlib);
        decoder.read(&compressed[..1], &mut [], 0).unwrap();
        assert!(decoder.bits.nbits > 0);
        decoder.reset(format);
        let empty = wrap(&[3, 0], &[], format);
        assert_eq!(decoder.read(&empty, &mut [], 0), Ok((empty.len(), 0)));
        assert!(decoder.is_done());

        decoder.reset(Format::Zlib);
        decoder.read(&compressed, &mut [0; 2], 0).unwrap();
        assert!(decoder.queued_output.is_some());
        decoder.reset(format);
        assert_eq!(decoder.read(&empty, &mut [], 0), Ok((empty.len(), 0)));
        assert!(decoder.is_done());
    }
}

#[test]
fn bounded_output_accepts_exact_size_and_rejects_excess() {
    for (raw, data) in fixtures() {
        let compressed = wrap(&raw, &data, Format::Zlib);
        match decompress_to_vec_bounded(&compressed, data.len()) {
            Ok(output) => assert_eq!(output, data),
            Err(_) => panic!("exact output limit rejected"),
        }
        if !data.is_empty() {
            assert!(matches!(
                decompress_to_vec_bounded(&compressed, data.len() - 1),
                Err(BoundedDecompressionError::OutputTooLarge { .. })
            ));
        }
    }
}
