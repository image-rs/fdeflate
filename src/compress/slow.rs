use std::io::{self, Write};

use super::{BitWriter, Symbol};

use super::bt_matchfinder::BTreeMatchFinder;

pub(super) struct SlowCompressor {
    match_finder: BTreeMatchFinder,

    min_match: u8,
    skip_ahead_shift: u8,
    search_depth: u16,
    nice_length: u16,
    max_lazy: u16,
}

impl SlowCompressor {
    pub fn new() -> Self {
        Self {
            match_finder: BTreeMatchFinder::new(3),

            min_match: 4,
            skip_ahead_shift: 9,
            search_depth: 64,
            nice_length: 258,
            max_lazy: 32,
        }
    }

    pub fn compress<W: Write>(&mut self, writer: &mut BitWriter<W>, data: &[u8]) -> io::Result<()> {
        let mut ip = 0;

        let mut length = 0;
        let mut distance = 0;
        let mut match_start = 0;

        while ip < data.len() {
            let mut symbols = Vec::new();
            let mut num_symbols = 0;

            let mut last_match = ip;
            'outer: while symbols.len() < 16384 && ip + 8 < data.len() {
                let current = u64::from_le_bytes(data[ip..][..8].try_into().unwrap());

                if length == 0 {
                    // if current == 0 {
                    //     while ip > last_match && data[ip - 1] == 0 {
                    //         ip -= 1;
                    //     }

                    //     if ip == 0 || data[ip - 1] != 0 {
                    //         ip += 1;
                    //     }

                    //     symbols.push(Symbol::LiteralRun {
                    //         start: last_match as u32,
                    //         end: ip as u32,
                    //     });
                    //     num_symbols += ip - last_match;

                    //     let mut run_length = 0;
                    //     while ip < data.len() && data[ip] == 0 && run_length < 258 {
                    //         run_length += 1;
                    //         ip += 1;
                    //     }

                    //     symbols.push(Symbol::Backref {
                    //         length: run_length as u16,
                    //         distance: 1,
                    //         dist_sym: 0,
                    //     });
                    //     num_symbols += 1;

                    //     last_match = ip;

                    //     length = 0;
                    //     continue;
                    // }

                    (length, distance, match_start) =
                        self.match_finder.get_and_insert(&data, ip, current, 4);
                }

                if length >= 3 {
                    if
                    /*match_start + length as usize > ip + 1
                    && length < self.max_lazy
                    &&*/
                    ip + length as usize + 9 <= data.len() {
                        ip += 1;
                        let (next_length, next_distance, next_match_start) = self
                            .match_finder
                            .get_and_insert(&data, ip, current >> 8, length + 1);
                        if next_length > 0 && match_start + 1 >= next_match_start {
                            assert!(next_length > length);
                            distance = next_distance;
                            length = next_length;
                            match_start = next_match_start;
                            continue;
                        }
                    }
                    assert!(last_match <= match_start);

                    symbols.push(Symbol::LiteralRun {
                        start: last_match as u32,
                        end: match_start as u32,
                    });
                    num_symbols += match_start - last_match;

                    symbols.push(Symbol::Backref {
                        length: length as u16,
                        distance,
                        dist_sym: super::distance_to_dist_sym(distance),
                    });
                    num_symbols += 1;

                    let match_end = match_start + length as usize;

                    if match_end + 8 < data.len() {
                        for j in (ip + 1)..match_end {
                            let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                            self.match_finder.insert(data, v, j);
                        }
                    }

                    ip = match_end;
                    last_match = match_end;

                    length = 0;
                    continue 'outer;
                }

                ip += 1;
            }
            if data.len() <= ip + 8 {
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
