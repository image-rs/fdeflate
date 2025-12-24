mod bitstream;
mod bitwriter;
mod matchfinder;
mod parse;
mod ultrafast;

use std::io::{self, Write};

use simd_adler32::Adler32;
pub use ultrafast::UltraFastCompressor;

use bitwriter::BitWriter;
use matchfinder::{HashChainMatchFinder, HashTableMatchFinder};
use parse::GreedyParser;

use crate::compress::parse::{LazyParser, RleParser};

const STORED_BLOCK_MAX_SIZE: usize = u16::MAX as usize;
const WINDOW_SIZE: usize = 32768;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
enum Flush {
    None,
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

/// Compressor that produces zlib or raw deflate compressed streams.
pub struct Compressor<W: Write> {
    /// State specific to the chosen compression level.
    inner: CompressorInner,

    /// Buffered input data.
    input: InputStream,

    /// Bit writer for writing output data.
    writer: BitWriter<W>,

    //// How much lookback to use for compression.
    window_size: usize,

    checksum: Option<Adler32>,
}
impl<W: Write> Compressor<W> {
    /// Create a new compressor.
    ///
    /// `level` is the compression level, where 0 is no compression values 1-9 are different levels
    /// of compression.
    ///
    /// `zlib` indicates whether to write a zlib header and checksum.
    pub fn new(mut writer: W, level: u8, zlib: bool) -> io::Result<Self> {
        if zlib {
            writer.write_all(&[0x78, 0x01])?; // zlib header
        }

        use CompressorInner::*;
        let inner = match level {
            0 => Uncompressed,
            1 => Fast(GreedyParser::new(5, HashTableMatchFinder::new())),
            2 => MediumFast(GreedyParser::new(6, HashChainMatchFinder::new(8, 16, 64))),
            3 => Medium(GreedyParser::new(6, HashChainMatchFinder::new(6, 16, 32))),
            4 => High(LazyParser::new(9, 12, HashChainMatchFinder::new(5, 16, 32))),
            5 => High(LazyParser::new(9, 16, HashChainMatchFinder::new(5, 32, 64))),
            6 => High(LazyParser::new(
                9,
                16,
                HashChainMatchFinder::new(4, 128, 128),
            )),
            7.. => High(LazyParser::new(
                12,
                128,
                HashChainMatchFinder::new(4, 512, 258),
            )),
        };

        Ok(Self {
            inner,
            writer: BitWriter::new(writer),
            input: InputStream {
                data: Vec::new(),
                base_index: 0,
                written: 0,
            },
            window_size: if level == 0 { 0 } else { WINDOW_SIZE },
            checksum: zlib.then(Adler32::new),
        })
    }

    /// Create a new compressor that only emits RLE matches.
    ///
    /// This mode is expected to be significantly faster even than level 1, but get somewhat worse
    /// compression ratios. It corresponds to the Z_RLE option in zlib.
    pub fn new_rle(mut writer: W, zlib: bool) -> io::Result<Self> {
        if zlib {
            writer.write_all(&[0x78, 0x01])?; // zlib header
        }

        Ok(Self {
            inner: CompressorInner::Rle(RleParser::new(5)),
            writer: BitWriter::new(writer),
            input: InputStream {
                data: Vec::new(),
                base_index: 0,
                written: 0,
            },
            window_size: 1,
            checksum: zlib.then(Adler32::new),
        })
    }

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> std::io::Result<()> {
        // Encoders use 32-bit indices in various places for performance. Limiting the input size
        // here simplifies things.
        const CHUNK_SIZE: usize = 1024 * 1024 * 1024;
        if data.len() > CHUNK_SIZE {
            for chunk in data.chunks(CHUNK_SIZE) {
                self.write_data(chunk)?;
            }
            return Ok(());
        }

        if let Some(ref mut checksum) = self.checksum {
            checksum.write(data);
        }

        // If we have no buffered input, attempt to compress directly from the input buffer.
        if self.input.data.is_empty() {
            let written = self.inner.compress(
                &mut self.writer,
                data,
                self.input.base_index,
                0,
                Flush::None,
            )?;

            let start = written.saturating_sub(self.window_size);
            self.input.data.extend_from_slice(&data[start..]);
            self.input.base_index += start as u32;
            self.input.written = written - start;
            return Ok(());
        }

        // If the indices used by the compressor would overflow, reset the base index.
        if u64::from(self.input.base_index) + data.len() as u64 > u64::from(u32::MAX) {
            match &mut self.inner {
                CompressorInner::Uncompressed => {}
                CompressorInner::Rle(rle) => rle.reset_indices(self.input.base_index),
                CompressorInner::Fast(fast) => fast.reset_indices(self.input.base_index),
                CompressorInner::MediumFast(medium) => medium.reset_indices(self.input.base_index),
                CompressorInner::Medium(medium_high) => {
                    medium_high.reset_indices(self.input.base_index)
                }
                CompressorInner::High(high) => high.reset_indices(self.input.base_index),
            }
            self.input.base_index = 0;
        }

        // Append the new data to the input buffer and compress it.
        self.input.data.extend_from_slice(data);
        let written = self.inner.compress(
            &mut self.writer,
            &self.input.data,
            self.input.base_index,
            self.input.written,
            Flush::None,
        )?;
        self.input.written += written;

        // Discard input data from before the start of the window, but avoid doing so too often.
        let discard = self.input.written.saturating_sub(self.window_size);
        if discard > 128 * 1024 {
            self.input.discard_bytes(discard);
        }

        Ok(())
    }

    /// Write the remainder of the stream and return the inner writer.
    pub fn finish(mut self) -> std::io::Result<W> {
        let written = self.inner.compress(
            &mut self.writer,
            &self.input.data,
            self.input.base_index,
            self.input.written,
            Flush::Finish,
        )?;

        self.input.data.clear();
        self.input.base_index += (self.input.written + written) as u32;
        self.input.written = 0;

        let writer = self.writer.flush()?;
        if let Some(checksum) = self.checksum.take() {
            let checksum_value = checksum.finish();
            writer.write_all(&checksum_value.to_be_bytes())?;
        }

        Ok(self.writer.take())
    }
}

enum CompressorInner {
    Uncompressed,
    Rle(RleParser),
    Fast(GreedyParser<HashTableMatchFinder>),
    MediumFast(GreedyParser<HashChainMatchFinder<true>>),
    Medium(GreedyParser<HashChainMatchFinder>),
    High(LazyParser<HashChainMatchFinder>),
}
impl CompressorInner {
    fn compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        input: &[u8],
        base_index: u32,
        start: usize,
        flush: Flush,
    ) -> std::io::Result<usize> {
        if flush == Flush::Finish && input.len() == start {
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

        if flush == Flush::Sync {
            writer.write_bits(0, 3)?;
            writer.flush()?.write_all(&[0, 0, 0xff, 0xff])?;
        }

        Ok(written)
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
