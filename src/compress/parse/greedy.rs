use std::io::{self, Write};

use crate::compress::{
    matchfinder::{Match, MatchFinder},
    parse::ParserInner,
    BitWriter, Flush,
};

pub(crate) struct GreedyParser<M> {
    inner: ParserInner<M>,
    m: Match,
}

impl<M: MatchFinder> GreedyParser<M> {
    pub fn new(skip_ahead_shift: u8, match_finder: M) -> Self {
        Self {
            inner: ParserInner::new(skip_ahead_shift, match_finder),
            m: Match::empty(),
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
        if !self.m.is_empty() {
            self.m.start -= delta;
        }

        let lookahead = if flush == Flush::None { 258 + 8 } else { 7 };
        let max_ip = data.len().saturating_sub(lookahead);

        loop {
            // Find a match if we don't yet have one.
            if self.m.is_empty() {
                self.m = self.inner.advance_to_match(data, base_index, max_ip);
                if self.m.is_empty() {
                    break;
                }
            }

            // Insert match finder entries for the current match.
            self.inner.advance(data, base_index, self.m.end());

            // Do a lookup at the position following the match. We'll need this even if we
            // accept the match, so it doesn't cost anything.
            let mut m2 = Match::empty();
            if self.inner.ip < max_ip {
                m2 = self.inner.get_match(data, base_index, true);
            } else if flush == Flush::None {
                break;
            }

            // Insert the current match, unless the next match starts too close to the current
            // one. Because we expand matches backwards, the next match might almost completely
            // overlap. If so, it'll probably be cheaper to emit an extra literal rather than an
            // extra backref.
            if m2.is_empty() || m2.start > self.m.start + 1 {
                self.inner.insert_match(base_index, &self.m);
                self.inner
                    .write_block_if_ready(writer, data, base_index, flush)?;

                // If the next match starts before the end of the current match, we need to
                // adjust the next match length and start position.
                if !m2.is_empty() && m2.start < self.inner.last_match {
                    assert!(m2.length >= 3);
                    m2.length -= (self.inner.last_match - m2.start) as u16;
                    m2.start = self.inner.last_match;
                    if m2.length < 4 {
                        m2 = Match::empty();
                    }
                }
            }

            // Advance to the next match (which might have a length of zero)
            self.m = m2;
        }

        self.inner
            .end_compress(writer, data, base_index, start, flush)
    }
}
