use std::io::{self, Write};

use super::{BitWriter, HashChainMatchFinder, Symbol};

pub(super) struct MediumCompressor {
    match_finder: HashChainMatchFinder,
    skip_ahead_shift: u8,
}

impl MediumCompressor {
    pub fn new(search_depth: u16, nice_length: u16, skip_ahead_shift: u8) -> Self {
        Self {
            match_finder: HashChainMatchFinder::new(search_depth, nice_length, 4),
            skip_ahead_shift,
        }
    }

    pub fn compress<W: Write>(&mut self, writer: &mut BitWriter<W>, data: &[u8]) -> io::Result<()> {
        let mut ip = 0; // Points at the next byte to hash/lookup for.
        let mut last_match = 0; //ip;

        while ip < data.len() {
            let mut length = 0;
            let mut distance = 0;
            let mut match_start = 0;

            let mut symbols = Vec::new();
            while symbols.len() < 16384 && ip + 8 <= data.len() {
                if length == 0 {
                    let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                    if current & 0xFF_FFFF_FFFF == 0 {
                        while ip > last_match && data[ip - 1] == 0 {
                            ip -= 1;
                        }

                        if ip == 0 || data[ip - 1] != 0 {
                            ip += 1;
                        }

                        symbols.push(Symbol::LiteralRun {
                            start: last_match as u32,
                            end: ip as u32,
                        });

                        let mut run_length = 0;
                        while ip < data.len() && data[ip] == 0 && run_length < 258 {
                            run_length += 1;
                            ip += 1;
                        }

                        symbols.push(Symbol::Backref {
                            length: run_length as u16,
                            distance: 1,
                            dist_sym: 0,
                        });

                        last_match = ip;

                        length = 0;
                        continue;
                    }

                    (length, distance, match_start) = self
                        .match_finder
                        .get_and_insert(&data, last_match, ip, current, 3);
                    ip += 1;
                }

                // If we haven't found a match in a while, start skipping ahead by emitting multiple
                // literals at once.
                if length < 3 {
                    ip += (ip - last_match) >> self.skip_ahead_shift;
                    continue;
                }

                assert!(last_match <= ip);
                assert!(last_match <= match_start,);
                let (mut next_length, mut next_distance, mut next_match_start) = (0, 0, 0);

                let match_end = match_start + length as usize;
                if match_end > ip {
                    // Insert match finder entries for the current match.
                    let insert_end = (match_end - 2).min(data.len() - 8);
                    let insert_start = ip.max(insert_end.saturating_sub(16));
                    for j in (insert_start..insert_end).step_by(3) {
                        let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                        self.match_finder.insert(v, j);
                        self.match_finder.insert(v >> 8, j + 1);
                        self.match_finder.insert(v >> 16, j + 2);
                    }
                    ip = match_end;

                    // Do a lookup at the position following the match. We'll need this even if we
                    // accept the match, so it doesn't cost anything.
                    if ip + 8 <= data.len() {
                        let next = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                        (next_length, next_distance, next_match_start) = self
                            .match_finder
                            .get_and_insert(&data, last_match, ip, next, 3);
                        ip += 1;
                    }
                }

                // Insert the current match, unless the next match starts too close to the current
                // one. Because we expand matches backwards, the next match might almost completely
                // overlap. If so, it'll probably be cheaper to emit an extra literal rather than an
                // extra backref.
                if next_length < 3 || next_match_start > match_start + 1 {
                    assert!(last_match <= match_start);
                    symbols.push(Symbol::LiteralRun {
                        start: last_match as u32,
                        end: match_start as u32,
                    });
                    symbols.push(Symbol::Backref {
                        length: length as u16,
                        distance,
                        dist_sym: super::distance_to_dist_sym(distance),
                    });
                    last_match = match_start + length as usize;

                    // If the next match starts before the end of the current match, we need to
                    // adjust the next match length and start position.
                    if next_length > 0 && next_match_start < last_match {
                        assert!(next_length >= 3);
                        next_length -= (last_match - next_match_start) as u16;
                        next_match_start = last_match;
                        if next_length < 4 {
                            next_length = 0;
                        }
                    }
                }

                // Advance to the next match (which might have a length of zero)
                length = next_length;
                match_start = next_match_start;
                distance = next_distance;
            }
            if data.len() < ip + 8 {
                symbols.push(Symbol::LiteralRun {
                    start: last_match as u32,
                    end: data.len() as u32,
                });
                ip = data.len();
            }
            super::write_block(writer, data, &symbols, ip == data.len())?;
        }

        Ok(())
    }
}
