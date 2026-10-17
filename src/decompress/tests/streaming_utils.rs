//! Bounded streaming driver shared by regression tests and the streaming fuzzer.

#[cfg(test)]
use crate as fdeflate;

use fdeflate::{DecompressionError, Decompressor, Format};

#[derive(Debug)]
pub struct StreamResult {
    pub output: Vec<u8>,
    pub consumed: usize,
    pub done: bool,
}

/// Decode with independently scheduled input chunks and writable output sizes. Zero-sized
/// chunks deliberately exercise backpressure. Unconsumed input is always offered again.
/// `stream_end`, when known, catches overconsumption on *every* call, not just at completion.
pub fn stream(
    input: &[u8],
    format: Format,
    output_limit: usize,
    input_chunks: &[usize],
    output_chunks: &[usize],
    stream_end: Option<usize>,
) -> Result<StreamResult, DecompressionError> {
    assert!(input_chunks.iter().any(|&n| n > 0));
    assert!(output_chunks.iter().any(|&n| n > 0));
    let mut decoder = Decompressor::new_with_format(format);
    let mut output = vec![0xff; output_limit];
    let (mut input_pos, mut output_pos, mut input_end, mut chunk_index) = (0usize, 0usize, 0, 0);
    // Each productive call advances at least one input/output byte. Allow intervening zero
    // chunks, plus a final call that recognizes completion without consuming another byte.
    let max_calls =
        (input.len() + output_limit + 2) * (input_chunks.len() + output_chunks.len() + 1);
    for call in 0..max_calls {
        if input_pos == input_end {
            input_end = input_pos
                .saturating_add(input_chunks[chunk_index % input_chunks.len()])
                .min(input.len());
            chunk_index += 1;
        }
        let output_end = output_pos
            .saturating_add(output_chunks[call % output_chunks.len()])
            .min(output_limit);
        let (consumed, produced) = decoder.read(
            &input[input_pos..input_end],
            &mut output[..output_end],
            output_pos,
        )?;
        assert!(consumed <= input_end - input_pos);
        assert!(produced <= output_end - output_pos);
        input_pos += consumed;
        output_pos += produced;
        if let Some(end) = stream_end {
            assert!(input_pos <= end, "consumed bytes after the first stream");
        }
        assert!(
            decoder.is_done() || input_pos == input_end || output_pos == output_end,
            "read violated its progress guarantee"
        );
        // A scheduled zero-sized buffer is temporary backpressure, not exhaustion of the
        // driver's input or output budget. Give the next scheduled chunk a chance to progress.
        let no_more_input =
            input_pos == input.len() && (output_pos < output_end || output_pos == output_limit);
        let no_more_output = output_pos == output_limit && input_pos < input_end;
        let stalled = consumed == 0 && produced == 0;
        if decoder.is_done() || (stalled && (no_more_input || no_more_output)) {
            if decoder.is_done() {
                if let Some(end) = stream_end {
                    assert_eq!(input_pos, end);
                }
                assert_eq!(
                    decoder.read(&input[input_pos..], &mut output, output_pos)?,
                    (0, 0)
                );
            }
            output.truncate(output_pos);
            return Ok(StreamResult {
                output,
                consumed: input_pos,
                done: decoder.is_done(),
            });
        }
    }
    panic!("streaming decoder exceeded its progress bound");
}
