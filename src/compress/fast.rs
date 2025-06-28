use std::io::{self, Write};

use super::{BitWriter, HashTableMatchFinder, Symbol};

pub(super) struct FastCompressor {
    match_finder: HashTableMatchFinder,
    skip_ahead_shift: u8,
}

impl FastCompressor {
    pub fn new() -> Self {
        Self {
            match_finder: HashTableMatchFinder::new(),
            skip_ahead_shift: 6,
        }
    }

    pub fn compress<W: Write>(&mut self, writer: &mut BitWriter<W>, data: &[u8]) -> io::Result<()> {
        let mut ip = 0;

        while ip < data.len() {
            let mut symbols = Vec::new();

            let mut last_match = ip;
            'outer: while symbols.len() < 16384 && ip + 8 <= data.len() {
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

                    continue;
                }

                let (length, distance, match_start) = self
                    .match_finder
                    .get_and_insert(&data, last_match, ip, current, 4);

                if length >= 3 {
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
                    let insert_start = (ip + 1).max(insert_end.saturating_sub(16));
                    for j in (insert_start..insert_end).step_by(3) {
                        let v = u64::from_le_bytes(data[j..][..8].try_into().unwrap());
                        self.match_finder.insert(v, j);
                    }

                    ip = match_end;
                    last_match = ip;

                    continue 'outer;
                }

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
