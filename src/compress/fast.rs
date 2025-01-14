use std::io::{self, Write};

use super::{BitWriter, HashChainMatchFinder, Symbol};

pub(super) struct FastCompressor {
    match_finder: HashChainMatchFinder,

    min_match: u8,
    skip_ahead_shift: u8,
    search_depth: u16,
    nice_length: u16,
}

impl FastCompressor {
    pub fn new() -> Self {
        Self {
            match_finder: HashChainMatchFinder::new(80, 16, 4),

            min_match: 4,
            skip_ahead_shift: 9,
            search_depth: 64,
            nice_length: 258,
        }
    }

    pub fn compress<W: Write>(&mut self, writer: &mut BitWriter<W>, data: &[u8]) -> io::Result<()> {
        let mut ip = 0;

        let mut length = 0;
        let mut distance = 0;
        let mut match_start = 0;

        while ip < data.len() {
            let mut symbols = Vec::new();

            let mut last_match = ip;
            'outer: while symbols.len() < 16384 && ip + 8 <= data.len() {
                let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());

                if length == 0 {
                    if current == 0 {
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
                        .get_and_insert(&data, last_match, ip, current, 4);
                }

                if length >= 3 {
                    // if match_start + length as usize > ip + 1
                    //     && length < max_lazy
                    //     && ip + length as usize + 9 <= data.len()
                    // {
                    //     ip += 1;
                    //     let (next_length, next_distance, next_match_start) =
                    //         matches.get_and_insert(&data, last_match, ip, current >> 8, length + 1);
                    //     if next_length > 0 && match_start + 1 >= next_match_start {
                    //         distance = next_distance;
                    //         length = next_length;
                    //         match_start = next_match_start;

                    //         if next_match_start > match_start {
                    //             continue;
                    //         }
                    //     }
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

                    let match_end = match_start + length as usize;
                    let insert_end = (match_end - 2).min(data.len() - 8);
                    let insert_start = ip + 1;//.max(insert_end.saturating_sub(16));
                    for j in (insert_start..insert_end).step_by(3) {
                        let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                        self.match_finder.insert(v, j);
                        self.match_finder.insert(v >> 8, j + 1);
                        self.match_finder.insert(v >> 16, j + 2);
                    }

                    // if insert_end >= insert_start + 3 {
                    //     let v = u64::from_le_bytes(data[insert_end - 3..][..8].try_into().unwrap());
                    //     matches.insert(v, insert_end - 3);
                    //     matches.insert(v >> 8, insert_end - 2);
                    //     matches.insert(v >> 16, insert_end - 1);
                    // }

                    ip = match_end;
                    last_match = ip;

                    length = 0;
                    continue 'outer;
                }

                // matches.insert(current >> 8, ip + 1);
                // matches.insert(current >> 16, ip + 2);

                // If we haven't found a match in a while, start skipping ahead by emitting multiple
                // literals at once.
                ip += 1 + ((ip - last_match) >> self.skip_ahead_shift);
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
