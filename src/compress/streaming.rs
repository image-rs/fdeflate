//! Bounded storage and backpressure for the block encoder.

use simd_adler32::Adler32;
use std::{
    fmt,
    io::{self, Write},
    mem::MaybeUninit,
};

use super::{
    bitwriter::BitState,
    parse::{LITERAL_LIMIT, SYMBOL_CAPACITY},
    BitWriter, CompressorInner, Flush, InputStream, RleParser, WINDOW_SIZE,
};

const BATCH_SIZE: usize = 128 * 1024;
const INPUT_CAPACITY: usize = WINDOW_SIZE + BATCH_SIZE;
const STORED_OUTPUT_CAPACITY: usize = BATCH_SIZE + 5 * BATCH_SIZE.div_ceil(65535) + 64;
// Fixed Huffman codes cost at most 9 bits per literal and 31 per backref.
// Pending tokens plus one input batch bound all output from one parser call:
// a new backref costs at most 9 bits per input byte it represents, too.
// At most nine symbol-limited blocks, two literal-limited blocks and the final
// block can finish in a batch. 64 bytes cover their framing, buffered bits,
// alignment, flush markers and the wrapper checksum. Stored/dynamic encodings
// are selected only when their complete bit cost is no greater than fixed.
const OUTPUT_CAPACITY: usize =
    (9 * (LITERAL_LIMIT + BATCH_SIZE) + 31 * SYMBOL_CAPACITY).div_ceil(8) + 64;

#[cfg(test)]
mod tests;

/// Requested action after all input supplied to a call has been accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CompressionFlush {
    /// Buffer input until a block is ready.
    None,
    /// Make all input available to a decoder without byte alignment.
    Partial,
    /// Make all input available and align the stream to a byte boundary.
    Sync,
    /// Perform a sync flush and discard compression history.
    Full,
    /// Finish the stream, including its checksum when using zlib.
    Finish,
}

/// The reason a streaming compression call stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionStatus {
    /// All supplied input was accepted; more input is needed to form a block.
    NeedsInput,
    /// Supply more output space and pass any unconsumed input again.
    NeedsOutput,
    /// The requested flush has been emitted and delivered.
    FlushComplete,
    /// The final block and optional checksum have been delivered.
    StreamEnd,
}

/// Incorrect use of a streaming compressor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionError {
    /// New input was supplied after finishing the stream.
    Finished,
}

impl fmt::Display for CompressionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("input supplied after finishing compression")
    }
}
impl std::error::Error for CompressionError {}

/// Input and output accounting for one compression call.
#[derive(Debug)]
pub struct CompressionProgress {
    /// Bytes accepted from the input slice.
    pub consumed: usize,
    /// Bytes initialized in the output slice.
    pub produced: usize,
    /// Completion, backpressure, or an error.
    pub status: Result<CompressionStatus, CompressionError>,
}

trait Output {
    fn len(&self) -> usize;
    fn put(&mut self, offset: usize, bytes: &[u8]);
}
impl Output for [u8] {
    fn len(&self) -> usize {
        self.len()
    }
    fn put(&mut self, offset: usize, bytes: &[u8]) {
        self[offset..offset + bytes.len()].copy_from_slice(bytes);
    }
}
impl Output for [MaybeUninit<u8>] {
    fn len(&self) -> usize {
        self.len()
    }
    fn put(&mut self, offset: usize, bytes: &[u8]) {
        for (destination, &byte) in self[offset..offset + bytes.len()].iter_mut().zip(bytes) {
            destination.write(byte);
        }
    }
}

// Select a sink once per batch. The direct path requires room for the worst
// case, avoiding a split/copy/overflow decision on every encoded machine word.
struct Destination<'a, O: Output + ?Sized> {
    output: &'a mut O,
    produced: &'a mut usize,
}
impl<O: Output + ?Sized> Write for Destination<'_, O> {
    #[inline]
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.output.put(*self.produced, bytes);
        *self.produced += bytes.len();
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct Pending<'a>(&'a mut Vec<u8>);
impl Write for Pending<'_> {
    #[inline]
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        assert!(self.0.len() + bytes.len() <= self.0.capacity());
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// A raw DEFLATE or zlib encoder with bounded working storage.
///
/// Calls continue across block boundaries until input or output is exhausted.
/// Pending output is drained before accepting more input. No caller buffer is
/// retained, and no allocation grows with the length of the stream.
pub struct StreamingCompressor {
    inner: CompressorInner,
    input: InputStream,
    bits: BitState,
    pending: Vec<u8>,
    pending_pos: usize,
    checksum: Option<Adler32>,
    window_size: usize,
    requested: Option<CompressionFlush>,
    completed: CompressionFlush,
    finished: bool,
}

impl StreamingCompressor {
    /// Creates a streaming compressor at compression level 0 through 9.
    pub fn new(level: u8, zlib: bool) -> Self {
        Self::with_parser(
            CompressorInner::new(level),
            if level == 0 { 0 } else { WINDOW_SIZE },
            zlib,
        )
    }

    /// Creates a streaming compressor restricted to run-length matches.
    pub fn new_rle(zlib: bool) -> Self {
        Self::with_parser(CompressorInner::Rle(RleParser::new(5)), 1, zlib)
    }

    fn with_parser(inner: CompressorInner, window_size: usize, zlib: bool) -> Self {
        let mut this = Self {
            inner,
            input: InputStream {
                data: Vec::with_capacity(INPUT_CAPACITY),
                base_index: 0,
                written: 0,
            },
            bits: BitState::default(),
            pending: Vec::with_capacity(if window_size == 0 {
                STORED_OUTPUT_CAPACITY
            } else {
                OUTPUT_CAPACITY
            }),
            pending_pos: 0,
            checksum: zlib.then(Adler32::new),
            window_size,
            requested: None,
            completed: CompressionFlush::None,
            finished: false,
        };
        if zlib {
            this.pending.extend_from_slice(&[0x78, 0x01]);
        }
        this
    }

    /// Resets the stream and wrapper format while reusing working allocations.
    pub fn reset(&mut self, zlib: bool) {
        self.inner.reset();
        self.input.data.clear();
        self.input.base_index = 0;
        self.input.written = 0;
        self.bits = BitState::default();
        self.pending.clear();
        self.pending_pos = 0;
        self.checksum = zlib.then(Adler32::new);
        if zlib {
            self.pending.extend_from_slice(&[0x78, 0x01]);
        }
        self.requested = None;
        self.completed = CompressionFlush::None;
        self.finished = false;
    }

    /// Processes input into initialized output storage.
    ///
    /// Pass unconsumed input again, retaining the flush argument. Once all input
    /// has been accepted, a flush remains pending until delivered, even if later
    /// calls use `None`. Repeating a completed flush without new input is a no-op.
    pub fn process(
        &mut self,
        input: &[u8],
        output: &mut [u8],
        flush: CompressionFlush,
    ) -> CompressionProgress {
        self.process_into(input, output, flush)
    }

    /// Processes input without initializing unused output capacity.
    pub fn process_uninit(
        &mut self,
        input: &[u8],
        output: &mut [MaybeUninit<u8>],
        flush: CompressionFlush,
    ) -> CompressionProgress {
        self.process_into(input, output, flush)
    }

    fn request(&mut self, flush: CompressionFlush) {
        if flush > self.completed && !self.finished {
            self.requested = Some(self.requested.unwrap_or(CompressionFlush::None).max(flush));
        }
    }

    fn process_into<O: Output + ?Sized>(
        &mut self,
        input: &[u8],
        output: &mut O,
        flush: CompressionFlush,
    ) -> CompressionProgress {
        let mut consumed = 0;
        let mut produced = 0;
        if self.finished && !input.is_empty() {
            return CompressionProgress {
                consumed,
                produced,
                status: Err(CompressionError::Finished),
            };
        }
        // In particular, preserve an empty-input flush even when output is full.
        if input.is_empty() {
            self.request(flush);
        }
        let status = loop {
            let n = (self.pending.len() - self.pending_pos).min(output.len() - produced);
            output.put(
                produced,
                &self.pending[self.pending_pos..self.pending_pos + n],
            );
            self.pending_pos += n;
            produced += n;
            if self.pending_pos != self.pending.len() {
                break CompressionStatus::NeedsOutput;
            }
            self.pending.clear();
            self.pending_pos = 0;
            if self.finished {
                break CompressionStatus::StreamEnd;
            }

            if let Some(requested) = self.requested {
                if produced == output.len() {
                    break CompressionStatus::NeedsOutput;
                }
                let action = match requested {
                    CompressionFlush::None => unreachable!(),
                    CompressionFlush::Partial => Flush::Partial,
                    CompressionFlush::Sync => Flush::Sync,
                    CompressionFlush::Full => Flush::Full,
                    CompressionFlush::Finish => Flush::Finish,
                };
                self.emit(output, &mut produced, action);
                self.completed = requested;
                self.requested = None;
                continue;
            }
            if consumed == input.len() {
                break if flush != CompressionFlush::None && self.completed >= flush {
                    CompressionStatus::FlushComplete
                } else {
                    CompressionStatus::NeedsInput
                };
            }
            if produced == output.len() {
                break CompressionStatus::NeedsOutput;
            }

            let buffered = self.input.data.len() - self.input.written;
            let n = (BATCH_SIZE - buffered).min(input.len() - consumed);
            if self.input.data.len() + n > INPUT_CAPACITY {
                // Rebase by a multiple of WINDOW_SIZE: hash-chain ring slots
                // must keep the same index modulo WINDOW_SIZE after rebasing.
                if self.input.base_index >= 1 << 30 {
                    self.rebase();
                }
                let discard = self.input.written.saturating_sub(self.window_size);
                self.input.discard_bytes(discard);
            }
            let fresh = &input[consumed..consumed + n];
            self.input.data.extend_from_slice(fresh);
            if let Some(checksum) = &mut self.checksum {
                checksum.write(fresh);
            }
            consumed += n;
            self.completed = CompressionFlush::None;
            if consumed == input.len() {
                self.request(flush);
            }
            if self.requested.is_none() && self.input.data.len() - self.input.written == BATCH_SIZE
            {
                self.emit(output, &mut produced, Flush::Batch);
            }
        };
        CompressionProgress {
            consumed,
            produced,
            status: Ok(status),
        }
    }

    fn rebase(&mut self) {
        let shift = self.input.base_index / WINDOW_SIZE as u32 * WINDOW_SIZE as u32;
        self.inner.reset_indices(shift);
        self.input.base_index -= shift;
    }

    fn emit<O: Output + ?Sized>(&mut self, output: &mut O, produced: &mut usize, flush: Flush) {
        debug_assert!(self.pending.is_empty());
        let raw = self.input.data.len() - self.input.written;
        let bound = if self.window_size == 0 {
            raw + 5 * raw.div_ceil(65535) + 64
        } else {
            (self.inner.pending_bits() + 9 * raw).div_ceil(8) + 64
        };
        debug_assert!(bound <= OUTPUT_CAPACITY);
        if output.len() - *produced >= bound {
            encode_batch(
                &mut self.inner,
                &mut self.input,
                &mut self.bits,
                self.checksum.as_ref(),
                Destination { output, produced },
                flush,
            );
        } else {
            encode_batch(
                &mut self.inner,
                &mut self.input,
                &mut self.bits,
                self.checksum.as_ref(),
                Pending(&mut self.pending),
                flush,
            );
        }
        self.finished = flush == Flush::Finish;
        if flush == Flush::Full {
            self.input.data.clear();
            self.input.written = 0;
            self.input.base_index = 0;
        }
    }
}

// Emission cannot be interrupted by output exhaustion: both sinks have room
// for the whole bounded batch before invoking the block encoder.
fn encode_batch<W: Write>(
    inner: &mut CompressorInner,
    input: &mut InputStream,
    bits: &mut BitState,
    checksum: Option<&Adler32>,
    destination: W,
    flush: Flush,
) {
    let mut writer = BitWriter::with_state(destination, *bits);
    let written = inner
        .compress(
            &mut writer,
            &input.data,
            input.base_index,
            input.written,
            flush,
        )
        .expect("the reserved block destination cannot fail");
    input.written += written;
    debug_assert_eq!(input.written, input.data.len());
    if flush == Flush::Finish {
        let destination = writer.flush().expect("infallible destination");
        if let Some(checksum) = checksum {
            destination
                .write_all(&checksum.finish().to_be_bytes())
                .expect("infallible destination");
        }
    } else if flush == Flush::Batch {
        // Release complete bytes of any blocks emitted during this batch.
        writer.partial_flush().expect("infallible destination");
    }
    *bits = writer.state();
}
