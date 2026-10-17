mod bitstream;
mod bitwriter;
mod matchfinder;
mod parse;
mod streaming;

pub use streaming::{
    CompressionError, CompressionFlush, CompressionProgress, CompressionStatus, StreamingCompressor,
};
mod ultrafast;

use std::io::{self, Write};

pub use ultrafast::UltraFastCompressor;

use bitwriter::BitWriter;
use matchfinder::{HashChainMatchFinder, HashTableMatchFinder};
use parse::GreedyParser;

use crate::compress::{
    matchfinder::HybridMatchFinder,
    parse::{LazyParser, RleParser},
};

const STORED_BLOCK_MAX_SIZE: usize = u16::MAX as usize;
const WINDOW_SIZE: usize = 32768;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
enum Flush {
    None,
    Batch,
    Partial,
    Full,
    Sync,
    Finish,
}

struct InputStream {
    data: Vec<u8>,
    base_index: u32,
    written: usize,
}
impl InputStream {
    fn discard_bytes(&mut self, n: usize) {
        self.data.copy_within(n.., 0);
        self.data.truncate(self.data.len() - n);
        self.base_index += n as u32;
        self.written -= n;
    }
}

/// The supported flush kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FlushKind {
    /// Finish compressing pending data and output an empty fixed block.
    Partial,
    /// Finish compressing pending data and output an empty stored block.
    Sync,
    /// Same as [`FlushKind::Sync`], but also resets the compression dictionary so subsequent output
    /// doesn't depend on previous output.
    Full,
}

/// Compressor that writes a bounded streaming encoder's output to a writer.
pub struct Compressor<W: Write> {
    inner: StreamingCompressor,
    writer: W,
    output: Box<[u8]>,
    failed: bool,
}
impl<W: Write> Compressor<W> {
    /// Creates a raw or zlib compressor at compression level 0 through 9.
    pub fn new(writer: W, level: u8, zlib: bool) -> io::Result<Self> {
        let mut compressor = Self {
            inner: StreamingCompressor::new(level, zlib),
            writer,
            output: vec![0; 32 * 1024].into_boxed_slice(),
            failed: false,
        };
        // Preserve the writer API's eager zlib header.
        compressor.run(&[], CompressionFlush::None)?;
        Ok(compressor)
    }

    /// Creates a compressor that only emits RLE matches.
    pub fn new_rle(writer: W, zlib: bool) -> io::Result<Self> {
        let mut compressor = Self {
            inner: StreamingCompressor::new_rle(zlib),
            writer,
            output: vec![0; 32 * 1024].into_boxed_slice(),
            failed: false,
        };
        compressor.run(&[], CompressionFlush::None)?;
        Ok(compressor)
    }

    fn run(&mut self, mut input: &[u8], flush: CompressionFlush) -> io::Result<()> {
        if self.failed {
            return Err(io::Error::other(
                "compressor cannot resume after a writer error",
            ));
        }
        let result = (|| loop {
            let progress = self.inner.process(input, &mut self.output, flush);
            input = &input[progress.consumed..];
            self.writer.write_all(&self.output[..progress.produced])?;
            let status = progress.status.map_err(io::Error::other)?;
            if status != CompressionStatus::NeedsOutput {
                return Ok(());
            }
        })();
        self.failed = result.is_err();
        result
    }

    /// Writes all input, emitting blocks as the bounded input buffer fills.
    /// After a writer error the compressor must be discarded.
    pub fn write_data(&mut self, data: &[u8]) -> io::Result<()> {
        self.run(data, CompressionFlush::None)
    }

    /// Flushes pending input and output with the specified flush kind.
    pub fn flush(&mut self, flush: FlushKind) -> io::Result<()> {
        self.run(
            &[],
            match flush {
                FlushKind::Partial => CompressionFlush::Partial,
                FlushKind::Sync => CompressionFlush::Sync,
                FlushKind::Full => CompressionFlush::Full,
            },
        )?;
        let result = self.writer.flush();
        self.failed = result.is_err();
        result
    }

    /// Gets a shared reference to the wrapped writer.
    pub fn get_ref(&self) -> &W {
        &self.writer
    }

    /// Gets a mutable reference to the wrapped writer.
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.writer
    }

    /// Finishes the stream and returns the writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.run(&[], CompressionFlush::Finish)?;
        Ok(self.writer)
    }
}

enum CompressorInner {
    Uncompressed,
    Rle(RleParser),
    Fast(GreedyParser<HashTableMatchFinder>),
    MediumFast(GreedyParser<HashChainMatchFinder<true>>),
    Medium(GreedyParser<HashChainMatchFinder>),
    High(LazyParser),
}
impl CompressorInner {
    fn new(level: u8) -> Self {
        match level {
            0 => Self::Uncompressed,
            1 => Self::Fast(GreedyParser::new(5, HashTableMatchFinder::new())),
            2 => Self::MediumFast(GreedyParser::new(6, HashChainMatchFinder::new(8, 16, 64))),
            3 => Self::Medium(GreedyParser::new(6, HashChainMatchFinder::new(6, 16, 32))),
            4 => Self::High(LazyParser::new(9, 12, HybridMatchFinder::new(5, 16, 32))),
            5 => Self::High(LazyParser::new(9, 16, HybridMatchFinder::new(5, 64, 64))),
            6 => Self::High(LazyParser::new(9, 16, HybridMatchFinder::new(4, 128, 128))),
            _ => Self::High(LazyParser::new(
                12,
                256,
                HybridMatchFinder::new(4, 256, 258),
            )),
        }
    }

    fn pending_bits(&self) -> usize {
        match self {
            Self::Uncompressed => 0,
            Self::Rle(p) => p.pending_bits(),
            Self::Fast(p) => p.pending_bits(),
            Self::MediumFast(p) => p.pending_bits(),
            Self::Medium(p) => p.pending_bits(),
            Self::High(p) => p.pending_bits(),
        }
    }

    fn reset(&mut self) {
        match self {
            Self::Uncompressed => {}
            Self::Rle(p) => p.reset(),
            Self::Fast(p) => p.reset(),
            Self::MediumFast(p) => p.reset(),
            Self::Medium(p) => p.reset(),
            Self::High(p) => p.reset(),
        }
    }

    fn compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        input: &[u8],
        base_index: u32,
        start: usize,
        flush: Flush,
    ) -> std::io::Result<usize> {
        if matches!(self, Self::Uncompressed) && flush == Flush::Finish && input.len() == start {
            writer.write_bits(3, 10)?;
            writer.flush()?;
            return Ok(0);
        }

        let written = match self {
            CompressorInner::Uncompressed => {
                let mut input = &input[start..];
                let mut written = 0;

                while input.len() > STORED_BLOCK_MAX_SIZE {
                    writer.write_bits(0, 3)?;
                    let writer = writer.flush()?;
                    writer.write_all(&[0xff, 0xff, 0, 0])?;
                    writer.write_all(&input[..STORED_BLOCK_MAX_SIZE])?;
                    input = &input[STORED_BLOCK_MAX_SIZE..];
                    written += STORED_BLOCK_MAX_SIZE;
                }

                if input.len() == STORED_BLOCK_MAX_SIZE || flush != Flush::None {
                    if flush == Flush::Finish {
                        writer.write_bits(1, 3)?;
                    } else {
                        writer.write_bits(0, 3)?;
                    }

                    let writer = writer.flush()?;
                    writer.write_all(&(input.len() as u16).to_le_bytes())?;
                    writer.write_all(&(!input.len() as u16).to_le_bytes())?;
                    writer.write_all(input)?;
                    written += input.len();
                }
                written
            }
            CompressorInner::Rle(rle) => rle.compress(writer, input, base_index, start, flush)?,
            CompressorInner::Fast(fast) => {
                fast.compress(writer, input, base_index, start, flush)?
            }
            CompressorInner::MediumFast(medium) => {
                medium.compress(writer, input, base_index, start, flush)?
            }
            CompressorInner::Medium(medium_high) => {
                medium_high.compress(writer, input, base_index, start, flush)?
            }
            CompressorInner::High(high) => {
                high.compress(writer, input, base_index, start, flush)?
            }
        };

        if flush == Flush::Partial {
            writer.write_bits(0b010, 10)?;
            writer.partial_flush()?;
        } else if flush == Flush::Full || flush == Flush::Sync {
            writer.write_bits(0, 3)?;
            let inner_writter = writer.flush()?;
            inner_writter.write_all(&[0, 0, 0xff, 0xff])?;
            inner_writter.flush()?;
        }

        Ok(written)
    }

    fn reset_indices(&mut self, old_base_index: u32) {
        match self {
            CompressorInner::Uncompressed => {}
            CompressorInner::Rle(rle) => rle.reset_indices(old_base_index),
            CompressorInner::Fast(fast) => fast.reset_indices(old_base_index),
            CompressorInner::MediumFast(medium) => medium.reset_indices(old_base_index),
            CompressorInner::Medium(medium_high) => medium_high.reset_indices(old_base_index),
            CompressorInner::High(high) => high.reset_indices(old_base_index),
        }
    }
}

/// Compresses the given data.
pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    compress_to_vec_with_level(input, 1)
}

/// Compresses the given data with a specific compression level.
pub fn compress_to_vec_with_level(input: &[u8], level: u8) -> Vec<u8> {
    let mut compressor = Compressor::new(Vec::with_capacity(input.len() / 4), level, true).unwrap();
    compressor.write_data(input).unwrap();
    compressor.finish().unwrap()
}

/// Compresses the given data using only RLE matches.
pub fn compress_to_vec_rle(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new_rle(Vec::with_capacity(input.len() / 4), true).unwrap();
    compressor.write_data(input).unwrap();
    compressor.finish().unwrap()
}

/// Compresses the given data using the ultra fast compression method.
pub fn compress_to_vec_ultra_fast(input: &[u8]) -> Vec<u8> {
    let mut compressor = UltraFastCompressor::new(Vec::with_capacity(input.len() / 4)).unwrap();
    compressor.write_data(input).unwrap();
    compressor.finish().unwrap()
}
