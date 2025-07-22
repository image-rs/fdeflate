mod bitwriter;
mod ultrafast;

use std::io::{self, Write};

use simd_adler32::Adler32;
pub use ultrafast::UltraFastCompressor;

use bitwriter::BitWriter;

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

        Ok(Self {
            inner: match level {
                0 => CompressorInner::Uncompressed,
                _ => CompressorInner::Uncompressed, // Placeholder for other compression levels
            },
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

    /// Write data to the compressor.
    pub fn write_data(&mut self, data: &[u8]) -> std::io::Result<()> {
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
}
impl CompressorInner {
    fn compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        input: &[u8],
        _base_index: u32,
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
                    writer.write_all(&input[start..][..STORED_BLOCK_MAX_SIZE])?;
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
