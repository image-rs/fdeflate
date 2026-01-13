use std::io::{self, Write};

use crate::compress::{
    matchfinder::{HybridMatchFinder, Match},
    parse::ParserInner,
    BitWriter, Flush,
};

pub(crate) struct LazyParser {
    inner: ParserInner<HybridMatchFinder>,
    max_lazy: u16,
    m0: Match,
    m1: Match,
}

impl LazyParser {
    pub fn new(skip_ahead_shift: u8, max_lazy: u16, match_finder: HybridMatchFinder) -> Self {
        Self {
            inner: ParserInner::new(skip_ahead_shift, match_finder),
            max_lazy,
            m1: Match::empty(),
            m0: Match::empty(),
        }
    }

    pub fn reset_indices(&mut self, old_base_index: u32) {
        self.inner.reset_indices(old_base_index);
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
        let delta = self.inner.start_compress(data, base_index, start);
        if !self.m0.is_empty() {
            self.m0.start -= delta;
        }
        if !self.m1.is_empty() {
            self.m1.start -= delta;
        }

        let lookahead = if flush == Flush::None { 258 + 8 } else { 7 };
        let max_ip = data.len().saturating_sub(lookahead);

        loop {
            if self.m1.is_empty() {
                self.m1 = self.inner.advance_to_match(data, base_index, max_ip);
                if self.m1.is_empty() {
                    break;
                }
            }

            assert!(self.inner.last_match <= self.inner.ip);
            assert!(self.inner.last_match <= self.m1.start);

            let mut m2 = Match::empty();
            if self.m1.length <= self.max_lazy {
                if self.inner.ip < max_ip {
                    assert!(self.inner.ip < self.m1.end());
                    let value = u64::from_le_bytes(data[self.inner.ip..][..8].try_into().unwrap());
                    m2 = self.inner.match_finder.get_and_insert_lazy(
                        data,
                        base_index,
                        self.inner.last_match,
                        self.inner.ip,
                        value,
                        self.m1.length + 1,
                    );
                    self.inner.ip += 1;
                    if m2.length <= self.m1.length {
                        m2 = Match::empty();
                    }
                } else if flush == Flush::None {
                    break;
                }
            }

            if m2.is_empty() {
                self.inner.advance(data, base_index, self.m1.end());

                if !self.m0.is_empty() && self.m0.start + 4 <= self.m1.start {
                    self.m0.length = self.m0.length.min((self.m1.start - self.m0.start) as u16);
                    self.inner.insert_match(base_index, self.m0);
                    self.m0 = Match::empty();
                }

                self.inner.insert_match(base_index, self.m1);
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

            self.inner
                .write_block_if_ready(writer, data, base_index, flush)?;
        }

        self.inner
            .end_compress(writer, data, base_index, start, flush)
    }
}
