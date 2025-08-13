use std::io::{self, Write};

use crate::compress::{
    bitstream::{self, Symbol},
    matchfinder::{Match, MatchFinder},
    BitWriter, Flush,
};

pub(crate) struct GreedyParser<M> {
    match_finder: M,
    skip_ahead_shift: u8,

    symbols: Vec<Symbol>,

    ip: usize,
    last_match: usize,
    m: Match,

    last_index: u32,
}

impl<M: MatchFinder> GreedyParser<M> {
    pub fn new(skip_ahead_shift: u8, match_finder: M) -> Self {
        Self {
            match_finder,
            skip_ahead_shift,
            symbols: Vec::new(),
            ip: 0,
            last_match: 0,
            m: Match::empty(),
            last_index: 0,
        }
    }

    fn rle_match(data: &[u8], last_match: usize, ip: usize) -> Match {
        let value = data[ip];

        let mut m = Match::new(4, 1, ip + 1);
        let min_start = 1.max(last_match).max(m.end().saturating_sub(258));

        while m.start > min_start && data[m.start - 2] == value {
            m.start -= 1;
            m.length += 1;
        }
        while m.length < 258 && data.get(m.end()) == Some(&value) {
            m.length += 1;
        }

        m
    }

    pub fn reset_indices(&mut self, old_base_index: u32) {
        self.last_index -= old_base_index;
        self.match_finder.reset_indices(old_base_index);
    }

    /// Compress the data using a greedy algorithm.
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
        if !self.m.is_empty() {
            self.m.start -= delta as usize;
        }
        self.last_index = base_index;

        let mut last_block_end = start;

        let max_ip = if flush == Flush::None {
            data.len().saturating_sub(258 + 8)
        } else {
            data.len().saturating_sub(7)
        };

        loop {
            if self.m.is_empty() {
                if self.ip >= max_ip {
                    break;
                }

                let current = u64::from_le_bytes(data[self.ip..][..8].try_into().unwrap());
                if current as u32 == (current >> 8) as u32 {
                    self.m = Self::rle_match(data, self.last_match, self.ip);
                    self.ip = self.m.end() - 3; // Skip inserting all the totally zero values into the hash table.
                } else {
                    self.m = self.match_finder.get_and_insert(
                        data,
                        base_index,
                        self.last_match,
                        self.ip,
                        current,
                    );
                    self.ip += 1;
                }

                if self.m.is_empty() {
                    // If we haven't found a match in a while, start skipping ahead by emitting
                    // multiple literals at once.
                    self.ip += (self.ip - self.last_match) >> self.skip_ahead_shift;
                    continue;
                }
            }

            assert!(self.last_match <= self.ip);
            assert!(self.last_match <= self.m.start);

            // Insert match finder entries for the current match.
            assert!(self.m.end() >= self.ip);
            for j in self.ip..self.m.end().min(data.len() - 8) {
                let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                self.match_finder.insert(v, base_index + j as u32);
            }
            self.ip = self.m.end();

            // Do a lookup at the position following the match. We'll need this even if we
            // accept the match, so it doesn't cost anything.
            let mut m2 = Match::empty();
            if self.ip < max_ip {
                let next = u64::from_le_bytes(data[self.ip..][..8].try_into().unwrap());
                if next as u32 == (next >> 8) as u32 {
                    m2 = Self::rle_match(data, self.last_match, self.ip);
                    // Skip inserting all the totally zero values into the hash table.
                    self.ip = m2.end() - 3;
                } else {
                    m2 = self
                        .match_finder
                        .get_and_insert(data, base_index, self.ip, self.ip, next);

                    while m2.length < 258
                        && m2.start > self.last_match
                        && m2.start > m2.distance as usize + 1
                        && data[m2.start - 1] == data[m2.start - m2.distance as usize - 1]
                    {
                        m2.length += 1;
                        m2.start -= 1;
                    }

                    self.ip += 1;
                }
            } else if flush == Flush::None {
                break;
            }

            // Insert the current match, unless the next match starts too close to the current
            // one. Because we expand matches backwards, the next match might almost completely
            // overlap. If so, it'll probably be cheaper to emit an extra literal rather than an
            // extra backref.
            if m2.is_empty() || m2.start > self.m.start + 1 {
                assert!(self.last_match <= self.m.start);
                if self.m.start > self.last_match {
                    self.symbols.push(Symbol::LiteralRun {
                        start: base_index + self.last_match as u32,
                        end: base_index + self.m.start as u32,
                    });
                }
                self.symbols.push(Symbol::Backref {
                    length: self.m.length,
                    distance: self.m.distance,
                    dist_sym: bitstream::distance_to_dist_sym(self.m.distance),
                });
                self.last_match = self.m.end();

                // If the next match starts before the end of the current match, we need to
                // adjust the next match length and start position.
                if m2.length > 0 && m2.start < self.last_match {
                    assert!(m2.length >= 3);
                    m2.length -= (self.last_match - m2.start) as u16;
                    m2.start = self.last_match;
                    if m2.length < 4 {
                        m2 = Match::empty();
                    }
                }
            }

            // Advance to the next match (which might have a length of zero)
            self.m = m2;

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
