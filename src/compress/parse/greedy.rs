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

    /// Bytes that have been encoded into symbols but not yet written to the output.
    pending_bytes: usize,
}

impl<M: MatchFinder> GreedyParser<M> {
    pub fn new(skip_ahead_shift: u8, match_finder: M) -> Self {
        Self {
            match_finder,
            skip_ahead_shift,
            symbols: Vec::new(),
            pending_bytes: 0,
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

        let mut ip = start + self.pending_bytes; // Points at the next byte to hash/lookup for.

        // If the last symbol was a literal run, remove it and set the last match position to be
        // before it.
        let mut last_match = ip;
        if let Some(Symbol::LiteralRun { start, end }) = self.symbols.last() {
            last_match -= (end - start) as usize;
            self.symbols.pop();
        }

        let mut m = Match::empty();
        let mut last_block_end = start;

        while ip + 8 <= data.len() {
            if m.is_empty() {
                let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                if current as u32 == (current >> 8) as u32 {
                    m = Self::rle_match(data, last_match, ip);
                    ip = m.end() - 3; // Skip inserting all the totally zero values into the hash table.
                } else {
                    m = self
                        .match_finder
                        .get_and_insert(&data, base_index, last_match, ip, current);
                    ip += 1;
                }

                if m.is_empty() {
                    // If we haven't found a match in a while, start skipping ahead by emitting
                    // multiple literals at once.
                    ip += (ip - last_match) >> self.skip_ahead_shift;
                    continue;
                }
            }

            assert!(last_match <= ip);
            assert!(last_match <= m.start);

            // Insert match finder entries for the current match.
            assert!(m.end() >= ip);
            for j in ip..m.end().min(data.len() - 8) {
                let v = u32::from_le_bytes(data[j..][..4].try_into().unwrap());
                self.match_finder.insert(v as u64, base_index + j as u32);
            }
            ip = m.end();

            // Do a lookup at the position following the match. We'll need this even if we
            // accept the match, so it doesn't cost anything.
            let mut m2 = Match::empty();
            if ip + 8 <= data.len() {
                let next = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                if next as u32 == (next >> 8) as u32 {
                    m2 = Self::rle_match(data, last_match, ip);
                    // Skip inserting all the totally zero values into the hash table.
                    ip = m2.end() - 3;
                } else {
                    m2 = self
                        .match_finder
                        .get_and_insert(&data, base_index, ip, ip, next);

                    while m2.length < 258
                        && m2.start > last_match
                        && m2.start > m2.distance as usize + 1
                        && data[m2.start - 1] == data[m2.start - m2.distance as usize - 1]
                    {
                        m2.length += 1;
                        m2.start -= 1;
                    }

                    ip += 1;
                }
            }

            // Insert the current match, unless the next match starts too close to the current
            // one. Because we expand matches backwards, the next match might almost completely
            // overlap. If so, it'll probably be cheaper to emit an extra literal rather than an
            // extra backref.
            if m2.is_empty() || m2.start > m.start + 1 {
                assert!(last_match <= m.start);
                if m.start > last_match {
                    self.symbols.push(Symbol::LiteralRun {
                        start: base_index + last_match as u32,
                        end: base_index + m.start as u32,
                    });
                }
                self.symbols.push(Symbol::Backref {
                    length: m.length as u16,
                    distance: m.distance,
                    dist_sym: bitstream::distance_to_dist_sym(m.distance),
                });
                last_match = m.end();

                // If the next match starts before the end of the current match, we need to
                // adjust the next match length and start position.
                if m2.length > 0 && m2.start < last_match {
                    assert!(m2.length >= 3);
                    m2.length -= (last_match - m2.start) as u16;
                    m2.start = last_match;
                    if m2.length < 4 {
                        m2 = Match::empty();
                    }
                }
            }

            // Advance to the next match (which might have a length of zero)
            m = m2;

            // Write the block if we have enough symbols.
            if self.symbols.len() >= 16384 {
                let last_block = flush == Flush::Finish && last_match == data.len();
                bitstream::write_block(writer, data, base_index, &self.symbols, last_block)?;
                self.symbols.clear();
                last_block_end = last_match;
            }
        }

        // The skip ahead logic can overshoot the end of the data.
        ip = ip.min(data.len());

        // Append a final literal run if there's remaining input.
        if last_match < data.len() {
            self.symbols.push(Symbol::LiteralRun {
                start: base_index + last_match as u32,
                end: base_index + data.len() as u32,
            });
            ip = data.len();
        }

        // If a flush was requested and there's remaining symbols, write another block.
        if flush != Flush::None && !self.symbols.is_empty() {
            assert_eq!(ip, data.len());

            let last_block = flush == Flush::Finish;
            bitstream::write_block(writer, data, base_index, &self.symbols, last_block)?;
            self.symbols.clear();
            last_block_end = ip;
        }

        self.pending_bytes = data.len() - last_block_end;
        Ok(last_block_end - start)
    }
}
