use std::io::{self, Write};

use crate::compress::{
    bitstream::{self, Symbol},
    matchfinder::{Match, MatchFinder},
    BitWriter, Flush,
};

pub(crate) struct LazyParser<M> {
    match_finder: M,
    skip_ahead_shift: u8,
    max_lazy: u16,

    symbols: Vec<Symbol>,

    ip: usize,
    last_match: usize,
    m0: Match,
    m1: Match,

    last_index: u32,
}

impl<M: MatchFinder> LazyParser<M> {
    pub fn new(skip_ahead_shift: u8, max_lazy: u16, match_finder: M) -> Self {
        Self {
            match_finder,
            skip_ahead_shift,
            max_lazy,
            symbols: Vec::new(),
            ip: 0,
            last_match: 0,
            m1: Match::empty(),
            m0: Match::empty(),
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

    fn find_match(&mut self, data: &[u8], base_index: u32, min_match: u16) -> Match {
        let value = u64::from_le_bytes(data[self.ip..][..8].try_into().unwrap());
        if value as u32 == (value >> 8) as u32 {
            let m = Self::rle_match(data, self.last_match, self.ip);
            self.ip = m.end() - 3; // Skip inserting all the totally zero values into the hash table.
            m
        } else {
            self.ip += 1;
            self.match_finder.get_and_insert(
                data,
                base_index,
                self.last_match,
                self.ip - 1,
                value,
                min_match,
            )
        }
    }

    fn accept_match(&mut self, m: Match) {
        assert!(self.last_match <= m.start);
        if m.start > self.last_match {
            self.symbols.push(Symbol::LiteralRun {
                start: self.last_index + self.last_match as u32,
                end: self.last_index + m.start as u32,
            });
        }
        self.symbols.push(Symbol::Backref {
            length: m.length,
            distance: m.distance,
            dist_sym: bitstream::distance_to_dist_sym(m.distance),
        });
        self.last_match = m.end();
    }

    fn advance(&mut self, data: &[u8], base_index: u32, new_ip: usize) {
        for j in self.ip..new_ip.min(data.len() - 8) {
            let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
            self.match_finder.insert(v, base_index + j as u32);
        }
        self.ip = new_ip;
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
        if !self.m1.is_empty() {
            self.m1.start -= delta as usize;
        }
        self.last_index = base_index;

        let mut last_block_end = start;

        let max_ip = if flush == Flush::None {
            data.len().saturating_sub(258 + 8)
        } else {
            data.len().saturating_sub(7)
        };

        loop {
            if self.m1.is_empty() {
                if self.ip >= max_ip {
                    break;
                }

                self.m1 = self.find_match(data, base_index, 4);
                if self.m1.is_empty() {
                    // If we haven't found a match in a while, start skipping ahead by emitting
                    // multiple literals at once.
                    self.ip += (self.ip - self.last_match) >> self.skip_ahead_shift;
                    continue;
                }
            }

            assert!(self.last_match <= self.ip);
            assert!(self.last_match <= self.m1.start);

            let mut m2 = Match::empty();
            if self.m1.length <= self.max_lazy {
                if self.ip < max_ip {
                    assert!(self.ip < self.m1.end());
                    m2 = self.find_match(data, base_index, self.m1.length + 1);
                    if m2.length <= self.m1.length {
                        m2 = Match::empty();
                    }
                } else if flush == Flush::None {
                    break;
                }
            }

            if m2.is_empty() {
                self.advance(data, base_index, self.m1.end());

                if !self.m0.is_empty() && self.m0.start + 4 <= self.m1.start {
                    self.m0.length = self.m0.length.min((self.m1.start - self.m0.start) as u16);
                    self.accept_match(self.m0);
                    self.m0 = Match::empty();
                }

                self.accept_match(self.m1);
                self.m0 = Match::empty();
                self.m1 = Match::empty();
                continue;
            } else if m2.start <= self.m1.start {
                self.m1 = m2;
                continue;
            } else {
                assert!(self.m1.length < m2.length);

                if self.m0.is_empty()
                    || self.m1.start < self.m0.start
                    || (self.m1.start == self.m0.start && self.m1.length > self.m0.length)
                {
                    self.m0 = self.m1;
                }
                self.m1 = m2;
            }

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
