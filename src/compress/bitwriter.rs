use std::io::{self, Write};

pub(super) struct BitWriter<W: Write> {
    buffer: u64,
    nbits: u8,
    writer: W,
}
impl<W: Write> BitWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            buffer: 0,
            nbits: 0,
            writer,
        }
    }

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

    pub fn flush(&mut self) -> io::Result<&mut W> {
        if self.nbits % 8 != 0 {
            self.write_bits(0, 8 - self.nbits % 8)?;
        }
        if self.nbits > 0 {
            self.writer
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])
                .unwrap();
            self.buffer = 0;
            self.nbits = 0;
        }
        Ok(&mut self.writer)
    }

    pub fn take(self) -> W {
        debug_assert_eq!(self.nbits, 0);
        self.writer
    }
}
