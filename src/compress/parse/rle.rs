use std::io::{self, Write};

use crate::compress::{
    bitstream::{self, Symbol},
    matchfinder::rle_match,
    BitWriter, Flush,
};

pub(crate) struct RleParser {
    skip_ahead_shift: u8,

    symbols: Vec<Symbol>,

    ip: usize,
    last_match: usize,

    last_index: u32,
}

impl RleParser {
    pub fn new(skip_ahead_shift: u8) -> Self {
        Self {
            skip_ahead_shift,
            symbols: Vec::new(),
            ip: 0,
            last_match: 0,
            last_index: 0,
        }
    }

    pub fn reset_indices(&mut self, old_base_index: u32) {
        self.last_index -= old_base_index;
    }

    pub fn compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        base_index: u32,
        start: usize,
        flush: Flush,
    ) -> io::Result<usize> {
        assert!(base_index as u64 + data.len() as u64 <= u32::MAX as u64);

        let delta = base_index - self.last_index;
        self.ip -= delta as usize;
        self.last_match -= delta as usize;
        self.last_index = base_index;

        let mut last_block_end = start;

        let max_ip = if flush == Flush::None {
            data.len().saturating_sub(258 + 7)
        } else {
            data.len().saturating_sub(7)
        };

        while self.ip < max_ip {
            let current = u64::from_le_bytes(data[self.ip..][..8].try_into().unwrap());
            if current as u32 != (current >> 8) as u32 {
                self.ip += 1 + ((self.ip - self.last_match) >> self.skip_ahead_shift);
                continue;
            }

            let m = rle_match(data, self.last_match, self.ip);
            self.ip = m.end();

            assert!(self.last_match <= self.ip);
            assert!(self.last_match <= m.start);
            if m.start > self.last_match {
                self.symbols.push(Symbol::LiteralRun {
                    start: base_index + self.last_match as u32,
                    end: base_index + m.start as u32,
                });
            }
            self.symbols.push(Symbol::Backref {
                length: m.length,
                distance: m.distance,
                dist_sym: bitstream::distance_to_dist_sym(m.distance),
            });
            self.last_match = m.end();

            // Write the block if we have enough symbols.
            if self.symbols.len() >= 16384 {
                let last_block = flush == Flush::Finish && self.last_match == data.len();
                bitstream::write_block(writer, data, base_index, &self.symbols, last_block)?;
                self.symbols.clear();
                last_block_end = self.last_match;
            }
        }

        // If a flush was requested and there's remaining symbols, write another block.
        if flush != Flush::None && (!self.symbols.is_empty() || self.last_match < data.len()) {
            // The skip ahead logic can overshoot the end of the data.
            self.ip = self.ip.min(data.len());

            // Append a final literal run if there's remaining input.
            if self.last_match < data.len() {
                self.symbols.push(Symbol::LiteralRun {
                    start: base_index + self.last_match as u32,
                    end: base_index + data.len() as u32,
                });
                self.ip = data.len();
                self.last_match = data.len();
            }
            assert_eq!(self.ip, data.len());

            let last_block = flush == Flush::Finish;
            bitstream::write_block(writer, data, base_index, &self.symbols, last_block)?;
            self.symbols.clear();
            last_block_end = self.ip;
        }

        Ok(last_block_end - start)
    }
}
