mod greedy;
mod lazy;
mod rle;

use std::io::{self, Write};

pub(crate) use greedy::*;
pub(crate) use lazy::*;
pub(crate) use rle::*;

pub(super) const SYMBOL_LIMIT: usize = 16384;
// A literal run takes one symbol, so its bytes need a separate storage limit.
// Keep the old batch size for literal-heavy input and allow long matches to
// span many batches without rebuilding Huffman trees.
pub(super) const LITERAL_LIMIT: usize = 128 * 1024;
// A lazy-parser iteration can append two literal/match pairs before checking.
pub(super) const SYMBOL_CAPACITY: usize = SYMBOL_LIMIT + 4;

use crate::compress::{
    bitstream::{self, Symbol},
    matchfinder::{rle_match, Match, MatchFinder},
    BitWriter, Flush,
};

struct ParserInner<M> {
    match_finder: M,
    skip_ahead_shift: u8,

    symbols: Vec<Symbol>,
    // Only literals retained from earlier batches are copied here. New runs
    // refer directly to the current input until save_literals is called.
    literals: Vec<u8>,
    literal_count: usize,
    saved_symbols: usize,
    block_len: usize,

    ip: usize,
    last_match: usize,
    last_index: u32,
    // Preserve match-search skipping across automatic block boundaries. This
    // score is capped independently of stream length and absolute indices.
    literal_streak: usize,
}
impl<M: MatchFinder> ParserInner<M> {
    fn new(skip_ahead_shift: u8, match_finder: M) -> Self {
        Self {
            match_finder,
            skip_ahead_shift,
            symbols: Vec::with_capacity(SYMBOL_CAPACITY),
            literals: Vec::with_capacity(LITERAL_LIMIT),
            literal_count: 0,
            saved_symbols: 0,
            block_len: 0,
            ip: 0,
            last_match: 0,
            last_index: 0,
            literal_streak: 0,
        }
    }

    fn reset_indices(&mut self, old_base_index: u32) {
        // Parser positions are relative to the input slice. Literal offsets
        // address owned block storage and are independent of stream indices.
        self.last_index -= old_base_index;
        self.match_finder.reset_indices(old_base_index);
    }

    fn reset(&mut self) {
        self.match_finder.clear();
        self.symbols.clear();
        self.literals.clear();
        self.literal_count = 0;
        self.saved_symbols = 0;
        self.block_len = 0;
        self.ip = 0;
        self.last_match = 0;
        self.last_index = 0;
        self.literal_streak = 0;
    }

    fn start_compress(&mut self, data: &[u8], base_index: u32) -> usize {
        assert!(base_index as u64 + data.len() as u64 <= u32::MAX as u64);

        let delta = base_index - self.last_index;
        self.ip -= delta as usize;
        self.last_match -= delta as usize;
        self.last_index = base_index;
        delta as usize
    }

    #[inline(always)]
    fn get_match(&mut self, data: &[u8], base_index: u32, fizzle: bool) -> Match {
        let current = u64::from_le_bytes(data[self.ip..][..8].try_into().unwrap());
        if current as u32 == (current >> 8) as u32 {
            let m = rle_match(data, self.last_match, self.ip);
            self.ip = m.end() - 3; // Skip inserting all the totally zero values into the hash table.
            m
        } else {
            let anchor = if fizzle { self.ip } else { self.last_match };
            let mut m = self
                .match_finder
                .get_and_insert(data, base_index, anchor, self.ip, current);
            if fizzle {
                while m.length < 258
                    && m.start > self.last_match
                    && m.start > m.distance as usize + 1
                    && data[m.start - 1] == data[m.start - m.distance as usize - 1]
                {
                    m.length += 1;
                    m.start -= 1;
                }
            }
            debug_assert!(m.is_empty() || self.last_match <= m.start);
            self.ip += 1;
            m
        }
    }

    #[inline(always)]
    fn advance_to_match(&mut self, data: &[u8], base_index: u32, max_ip: usize) -> Match {
        while self.ip < max_ip {
            let m = self.get_match(data, base_index, false);
            if !m.is_empty() {
                return m;
            }

            // If we haven't found a match in a while, start skipping ahead by emitting
            // multiple literals at once.
            self.ip += (self.ip - self.last_match + self.literal_streak) >> self.skip_ahead_shift;
        }

        Match::empty()
    }

    /// Insert match finder entries for the given range.
    #[inline(always)]
    fn advance(&mut self, data: &[u8], base_index: u32, end: usize) {
        assert!(self.last_match <= self.ip);

        for j in self.ip..end.min(data.len() - 8) {
            let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
            self.match_finder.insert(v, base_index + j as u32);
        }
        self.ip = self.ip.max(end);
    }

    #[inline(always)]
    fn insert_match<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        m: Match,
    ) -> io::Result<()> {
        assert!(self.last_match <= m.start);
        if m.start > self.last_match {
            self.append_literals(writer, data, m.start)?;
        }
        debug_assert!(self.symbols.len() < SYMBOL_CAPACITY);
        self.symbols.push(Symbol::Backref {
            length: m.length,
            distance: m.distance,
            dist_sym: bitstream::distance_to_dist_sym(m.distance),
        });
        self.block_len += m.length as usize;
        self.last_match = m.end();
        self.literal_streak = 0;
        Ok(())
    }

    // Keep literal bookkeeping and block emission out of the match-only hot path.
    fn append_literals<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        end: usize,
    ) -> io::Result<()> {
        while self.last_match < end {
            if self.literal_count == LITERAL_LIMIT || self.symbols.len() >= SYMBOL_LIMIT {
                self.write_pending_block(writer, data, false)?;
            }
            let n = (end - self.last_match).min(LITERAL_LIMIT - self.literal_count);
            // Offsets address saved literals followed by the current input.
            // Copy only if the block outlives this input batch.
            let start = self.literals.len() + self.last_match;
            self.symbols.push(Symbol::LiteralRun {
                start: start as u32,
                end: (start + n) as u32,
            });
            self.literal_count += n;
            self.block_len += n;
            self.last_match += n;
        }
        Ok(())
    }

    #[inline(always)]
    fn write_block_if_ready<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
    ) -> io::Result<()> {
        if self.symbols.len() >= SYMBOL_LIMIT {
            self.write_pending_block(writer, data, false)?;
        }
        Ok(())
    }

    fn write_pending_block<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        eof: bool,
    ) -> io::Result<()> {
        // The block may span many input-window slides. Stored encoding is only
        // available while its entire raw representation remains in the window.
        let raw = self
            .last_match
            .checked_sub(self.block_len)
            .map(|start| &data[start..self.last_match]);
        bitstream::write_block(writer, &self.literals, data, raw, &self.symbols, eof)?;
        self.symbols.clear();
        self.literals.clear();
        self.literal_count = 0;
        self.saved_symbols = 0;
        self.block_len = 0;
        Ok(())
    }

    fn save_literals(&mut self, data: &[u8]) {
        let saved = self.literals.len();
        for symbol in &mut self.symbols[self.saved_symbols..] {
            if let Symbol::LiteralRun { start, end } = symbol {
                let bytes = &data[*start as usize - saved..*end as usize - saved];
                *start = self.literals.len() as u32;
                self.literals.extend_from_slice(bytes);
                *end = self.literals.len() as u32;
            }
        }
        self.saved_symbols = self.symbols.len();
        debug_assert_eq!(self.literals.len(), self.literal_count);
        debug_assert!(self.literal_count <= LITERAL_LIMIT);
        debug_assert_eq!(self.literals.capacity(), LITERAL_LIMIT);
        debug_assert_eq!(self.symbols.capacity(), SYMBOL_CAPACITY);
    }

    fn pending_bits(&self) -> usize {
        9 * self.literal_count + 31 * self.symbols.len()
    }

    fn end_compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        start: usize,
        flush: Flush,
    ) -> io::Result<usize> {
        debug_assert!(flush != Flush::None);
        if self.last_match < data.len() {
            self.literal_streak =
                (self.literal_streak + data.len() - self.last_match).min(1024 * 1024);
            self.append_literals(writer, data, data.len())?;
        }
        self.ip = data.len();
        // A parsing batch can finish without closing the Huffman block.
        if flush != Flush::Batch
            || self.literal_count == LITERAL_LIMIT
            || self.symbols.len() >= SYMBOL_LIMIT
        {
            if !self.symbols.is_empty() {
                self.write_pending_block(writer, data, flush == Flush::Finish)?;
            } else if flush == Flush::Finish {
                writer.write_bits(3, 10)?;
            }
        }
        if flush == Flush::Batch {
            self.save_literals(data);
        } else {
            self.literal_streak = 0;
        }
        if flush == Flush::Full {
            self.match_finder.clear();
            self.ip = 0;
            self.last_match = 0;
            self.last_index = 0;
            assert!(self.symbols.is_empty());
        }
        Ok(data.len() - start)
    }
}
