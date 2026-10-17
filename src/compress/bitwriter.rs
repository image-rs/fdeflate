use std::io::{self, Write};

#[derive(Clone, Copy, Default)]
pub(super) struct BitState {
    buffer: u64,
    nbits: u8,
}

pub(super) struct BitWriter<W: Write> {
    buffer: u64,
    nbits: u8,
    writer: W,
}
impl<W: Write> BitWriter<W> {
    pub fn with_state(writer: W, state: BitState) -> Self {
        Self {
            buffer: state.buffer,
            nbits: state.nbits,
            writer,
        }
    }

    pub fn state(&self) -> BitState {
        BitState {
            buffer: self.buffer,
            nbits: self.nbits,
        }
    }

    pub fn buffered_bits(&self) -> u8 {
        self.nbits
    }
    #[inline(always)]
    pub fn write_bits(&mut self, bits: u64, nbits: u8) -> io::Result<()> {
        debug_assert!(nbits <= 64);

        self.buffer |= bits << self.nbits;
        self.nbits += nbits;

        if self.nbits >= 64 {
            self.writer.write_all(&self.buffer.to_le_bytes())?;
            self.nbits -= 64;
            self.buffer = bits.checked_shr((nbits - self.nbits) as u32).unwrap_or(0);
        }
        debug_assert!(self.nbits < 64);
        Ok(())
    }

    /// Flush buffered full bytes to the underlying writer, and then flush the
    /// underlying writer. Any partial bytes in the buffer are left there.
    pub fn partial_flush(&mut self) -> io::Result<()> {
        if self.nbits >= 8 {
            self.writer
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])?;
            self.buffer >>= (self.nbits / 8) * 8;
            self.nbits %= 8;
        }

        self.writer.flush()?;

        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<&mut W> {
        if !self.nbits.is_multiple_of(8) {
            self.write_bits(0, 8 - self.nbits % 8)?;
        }
        if self.nbits > 0 {
            self.writer
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])?;
            self.buffer = 0;
            self.nbits = 0;
        }
        Ok(&mut self.writer)
    }
}
