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
            let mut length = 0u16;
            let mut distance = 0;
            let mut match_start = 0;

            let mut symbols = Vec::new();
            while symbols.len() < 16384 && ip + 8 <= data.len() {
                if length == 0 {
                    let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                    // if current & 0xFF_FFFF_FFFF == 0 {
                    //     length = 4;
                    //     match_start = ip + 1;
                    //     distance = 1;

                    //     let min_start = 1.max(last_match).max(match_start.saturating_sub(258 - 4));

                    //     while match_start > min_start && data[match_start - 2] == 0 {
                    //         match_start -= 1;
                    //         length += 1;
                    //     }
                    //     while length < 258
                    //         && match_start + (length as usize) < data.len()
                    //         && data[match_start + length as usize] == 0
                    //     {
                    //         length += 1;
                    //     }

                    //     // Skip inserting all the totally zero values into the hash table.
                    //     ip = match_start + length as usize - 3;
                    // } else {
                        (length, distance, match_start) = self.match_finder.get_and_insert(
                            &data,
                            last_match,
                            ip,
                            current as u32,
                            3,
                        );
                        ip += 1;
                    // }
                }

                if length < 3 {
                    // If we haven't found a match in a while, start skipping ahead by emitting
                    // multiple literals at once.
                    ip += (ip - last_match) >> self.skip_ahead_shift;
                    continue;
                }

                assert!(last_match <= ip);
                assert!(last_match <= match_start,);
                let (mut next_length, mut next_distance, mut next_match_start) = (0, 0, 0);

                let match_end = match_start + length as usize;
                if match_end >= ip {
                    // // Insert match finder entries for the current match.
                    // let insert_end = (match_end - 3).min(data.len() - 8);
                    // let insert_start = ip.max(insert_end.saturating_sub(16));
                    // for j in (insert_start..insert_end).step_by(4) {
                    //     let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                    //     self.match_finder.insert(v, j);
                    //     self.match_finder.insert(v >> 8, j + 1);
                    //     self.match_finder.insert(v >> 16, j + 2);
                    //     self.match_finder.insert(v >> 24, j + 3);
                    // }
                    for j in ip..match_end.min(data.len() - 8) {
                        let v = u32::from_le_bytes(data[j..][..4].try_into().unwrap());
                        self.match_finder.insert(v as u64, j);
                    }

                    ip = match_end;

                    // innumerable::event!("current-delta", ip as i32 - match_start as i32);

                    // Do a lookup at the position following the match. We'll need this even if we
                    // accept the match, so it doesn't cost anything.
                    if ip + 8 <= data.len() {
                        let next = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());
                        // if next & 0xFF_FFFF_FFFF == 0 {
                        //     next_length = 4;
                        //     next_match_start = ip + 1;
                        //     next_distance = 1;

                        //     let min_start =
                        //         1.max(last_match).max(next_match_start.saturating_sub(258 - 4));

                        //     while next_match_start > min_start && data[next_match_start - 2] == 0 {
                        //         next_match_start -= 1;
                        //         next_length += 1;
                        //     }
                        //     while next_length < 258
                        //         && next_match_start + (next_length as usize) < data.len()
                        //         && data[next_match_start + next_length as usize] == 0
                        //     {
                        //         next_length += 1;
                        //     }

                        //     // Skip inserting all the totally zero values into the hash table.
                        //     ip = next_match_start + next_length as usize - 3;
                        // } else {
                            (next_length, next_distance, next_match_start) = self
                                .match_finder
                                .get_and_insert(&data, last_match, ip, next as u32, 3);

                            // innumerable::event!("x-delta", next_match_start as i32 - ip as i32);

                            ip += 1;
                        // }
                    }
                }

                // if next_length >= 3 {
                //     // innumerable::event!("next-length", next_length);
                //     innumerable::event!("next-delta", next_match_start as i32 - match_start as i32);
                // }

                // Insert the current match, unless the next match starts too close to the current
                // one. Because we expand matches backwards, the next match might almost completely
                // overlap. If so, it'll probably be cheaper to emit an extra literal rather than an
                // extra backref.
                if next_length < 3 || next_match_start > match_start + 1 {
                    // if next_length < 3 && next_match_start > match_start + 1 {
                    //     innumerable::event!("match", 0);
                    // } else if next_length < 3 {
                    //     innumerable::event!("match", 1);
                    // } else {
                    //     innumerable::event!("match", 2);
                    // }

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
                    //     innumerable::event!("fizzle", 0);
                    // } else if next_length >= 3 {
                    //     innumerable::event!("fizzle", 1);
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
