use crate::decompress::{EXCEPTIONAL_ENTRY, LITERAL_ENTRY, SECONDARY_TABLE_ENTRY};

/// Return the next code, or if the codeword is already all ones (which is the final code), return
/// the same code again.
fn next_codeword(mut codeword: u16, table_size: u16) -> u16 {
    if codeword == table_size - 1 {
        return codeword;
    }

    let adv = (u16::BITS - 1) - (codeword ^ (table_size - 1)).leading_zeros();
    let bit = 1 << adv;
    codeword &= bit - 1;
    codeword |= bit;
    codeword
}

#[allow(clippy::needless_range_loop)]
pub fn build_table(
    lengths: &[u8],
    entries: &[u32],
    codes: &mut [u16],
    primary_table: &mut [u32],
    secondary_table: &mut Vec<u16>,
    is_distance_table: bool,
    double_literal: bool,
) -> bool {
    // Count the number of symbols with each code length.
    let mut histogram = [0; 16];
    for &length in lengths {
        histogram[length as usize] += 1;
    }

    // Determine the maximum code length.
    let mut max_length = 15;
    while max_length > 1 && histogram[max_length] == 0 {
        max_length -= 1;
    }

    // Handle zero and one symbol huffman codes (which are only allowed for distance codes).
    if is_distance_table {
        if max_length == 0 {
            primary_table.fill(0);
            secondary_table.clear();
            return true;
        } else if max_length == 1 && histogram[1] == 1 {
            let symbol = lengths.iter().position(|&l| l == 1).unwrap();
            codes[symbol] = 0;
            let entry = entries
                .get(symbol)
                .cloned()
                .unwrap_or((symbol as u32) << 16)
                | 1;
            for chunk in primary_table.chunks_mut(2) {
                chunk[0] = entry;
                chunk[1] = 0;
            }
            return true;
        }
    }

    // Sort symbols by code length. Given the histogram, we can determine the starting offset
    // for each code length.
    let mut offsets = [0; 16];
    let mut codespace_used = 0;
    offsets[1] = histogram[0];
    for i in 1..max_length {
        offsets[i + 1] = offsets[i] + histogram[i];
        codespace_used = (codespace_used << 1) + histogram[i];
    }
    codespace_used = (codespace_used << 1) + histogram[max_length];

    // Check that the provided lengths form a valid Huffman tree.
    if codespace_used != (1 << max_length) {
        return false;
    }

    // Sort the symbols by code length.
    let mut next_index = offsets;
    let mut sorted_symbols = [0; 288];
    for symbol in 0..lengths.len() {
        let length = lengths[symbol];
        sorted_symbols[next_index[length as usize]] = symbol;
        next_index[length as usize] += 1;
    }

    let mut codeword = 0u16;
    let mut i = histogram[0];

    // Populate the primary decoding table
    let primary_table_bits = primary_table.len().ilog2() as usize;
    let primary_table_mask = (1 << primary_table_bits) - 1;
    for length in 1..=primary_table_bits {
        let current_table_end = 1 << length;

        // Loop over all symbols with the current code length and set their table entries.
        for _ in 0..histogram[length] {
            let symbol = sorted_symbols[i];
            i += 1;

            primary_table[codeword as usize] = entries
                .get(symbol)
                .cloned()
                .unwrap_or((symbol as u32) << 16)
                | length as u32;

            codes[symbol] = codeword;
            codeword = next_codeword(codeword, current_table_end as u16);
        }

        if double_literal {
            for len1 in 1..length {
                let len2 = length - len1;
                for sym1_index in offsets[len1]..next_index[len1] {
                    for sym2_index in offsets[len2]..next_index[len2] {
                        let sym1 = sorted_symbols[sym1_index];
                        let sym2 = sorted_symbols[sym2_index];
                        if sym1 < 256 && sym2 < 256 {
                            let codeword1 = codes[sym1];
                            let codeword2 = codes[sym2];
                            let codeword = codeword1 | (codeword2 << len1);
                            let entry = (sym1 as u32) << 16
                                | (sym2 as u32) << 24
                                | LITERAL_ENTRY
                                | (2 << 8);
                            primary_table[codeword as usize] = entry | (length as u32);
                        }
                    }
                }
            }
        }

        // If we aren't at the maximum table size, double the size of the table.
        if length < primary_table_bits {
            primary_table.copy_within(0..current_table_end, current_table_end);
        }
    }

    // Populate the secondary decoding table.
    secondary_table.clear();
    if max_length > primary_table_bits {
        let mut subtable_start = 0;
        let mut subtable_prefix = !0;
        for length in (primary_table_bits + 1)..=max_length {
            let subtable_size = 1 << (length - primary_table_bits);
            let overflow_bits_mask = subtable_size as u32 - 1;
            for _ in 0..histogram[length] {
                // If the codeword's prefix doesn't match the current subtable, create a new
                // subtable.
                if codeword & primary_table_mask != subtable_prefix {
                    subtable_prefix = codeword & primary_table_mask;
                    subtable_start = secondary_table.len();
                    primary_table[subtable_prefix as usize] = ((subtable_start as u32) << 16)
                        | EXCEPTIONAL_ENTRY
                        | SECONDARY_TABLE_ENTRY
                        | overflow_bits_mask;
                    secondary_table.resize(subtable_start + subtable_size, 0);
                }

                // Lookup the symbol.
                let symbol = sorted_symbols[i];
                i += 1;

                // Insert the symbol into the secondary table and advance to the next codeword.
                codes[symbol] = codeword;
                secondary_table[subtable_start + (codeword >> primary_table_bits) as usize] =
                    ((symbol as u16) << 4) | (length as u16);
                codeword = next_codeword(codeword, 1 << length);
            }

            // If there are more codes with the same subtable prefix, extend the subtable.
            if length < max_length && codeword & primary_table_mask == subtable_prefix {
                secondary_table.extend_from_within(subtable_start..);
                let subtable_size = secondary_table.len() - subtable_start;
                let overflow_bits_mask = subtable_size as u32 - 1;
                primary_table[subtable_prefix as usize] = ((subtable_start as u32) << 16)
                    | EXCEPTIONAL_ENTRY
                    | SECONDARY_TABLE_ENTRY
                    | overflow_bits_mask;
            }
        }
    }

    true
}

#[cfg(test)]
mod test {
    use super::{LITERAL_ENTRY, SECONDARY_TABLE_ENTRY};
    use crate::tables::LITLEN_TABLE_ENTRIES;

    fn validate_tables(
        primary_table_bits: usize,
        lengths: &[u8],
        primary_table: &[u32],
        secondary_table: &[u16],
    ) {
        let expecting_only_double_literals =
            (*lengths.iter().max().unwrap() as usize) * 2 <= primary_table_bits;
        for (i, entry) in primary_table.into_iter().enumerate() {
            if 0 != entry & LITERAL_ENTRY {
                // Expected format: aaaaaaaa_bbbbbbbb_100000yy_0000xxxx
                match entry >> 8 & 0x7f {
                    1 => {
                        if expecting_only_double_literals {
                            panic!(
                                "Unexpected single literal: index={i} ({i:b}); entry=0b{entry:b}"
                            );
                        }
                    }
                    2 => (),
                    other => panic!("Unexpected output_advance_bytes={other}: index={i} ({i:b})"),
                }

                let input_bits = entry & 0xff;
                if input_bits == 0 {
                    panic!("input_advance_bits unexpectedly equal to 0");
                } else if input_bits > 15 {
                    panic!("Unexpectedly big input_advance_bits: {}", input_bits);
                }

                let symbol_mask = (1 << lengths.len().min(256).ilog2() + 1) - 1;
                let s1 = entry >> 16 & 0xff;
                if 0 != s1 & !symbol_mask {
                    panic!("Unexpectedly big symbol: {}", s1);
                }
                let s2 = entry >> 24 & 0xff;
                if 0 != s2 & !symbol_mask {
                    panic!("Unexpectedly big symbol: {}", s2);
                }
            } else if 0 != entry & SECONDARY_TABLE_ENTRY {
                // Expected format: 0000xxxx_xxxxxxxx_01100000_mmmmmmmm
                let overflow_bits_mask = (entry & 0xff) as usize;
                let overflow_bits = overflow_bits_mask.trailing_ones() as usize;
                if overflow_bits == 0 {
                    panic!("Unexpectedly missing mask: index={i} ({i:b}), entry={entry:b}");
                }
                if overflow_bits + primary_table_bits > 15 {
                    // Section 3.2.7 of https://www.ietf.org/rfc/rfc1951.txt implies
                    // that codeword lengths are at most 15.
                    panic!("Unexpectedly long symbol: index={i} ({i:b}), entry={entry:b}");
                }
                let index2_base = (entry >> 16) as usize;
                assert!(index2_base + overflow_bits_mask <= secondary_table.len());
            } else {
                // TODO: Provide test coverage/support for EOF symbol (257th symbol - 256)
                // and distance codes (even bigger symbols).
                assert!(lengths.len() > 256);
            }
        }
    }

    #[derive(Debug, Eq, PartialEq)]
    enum LitlenResult {
        SingleLiteral { symbol: u8, input_bits: usize },
        DoubleLiteral { s1: u8, s2: u8, input_bits: usize },
        SecondaryTableLiteral { symbol: u16, input_bits: usize },
    }

    struct LitlenTables {
        primary_table_bits: usize,
        primary_table_mask: u64,
        primary_table: Vec<u32>,
        secondary_table: Vec<u16>,
    }

    impl LitlenTables {
        fn new(primary_table_bits: usize, lengths: &[u8]) -> Option<Self> {
            let primary_table_size = 1 << primary_table_bits;
            let primary_table_mask = (primary_table_size - 1).try_into().unwrap();
            let mut primary_table = vec![0; primary_table_size];
            let mut secondary_table = Vec::new();
            let mut codes = [0; 288];

            const IS_DISTANCE_TABLE: bool = false;
            const DOUBLE_LITERAL: bool = true;

            let success = super::build_table(
                lengths,
                &LITLEN_TABLE_ENTRIES,
                &mut codes,
                &mut primary_table,
                &mut secondary_table,
                IS_DISTANCE_TABLE,
                DOUBLE_LITERAL,
            );

            if success {
                validate_tables(
                    primary_table_bits,
                    lengths,
                    &primary_table,
                    &secondary_table,
                );
                Some(Self {
                    primary_table_bits,
                    primary_table_mask,
                    primary_table,
                    secondary_table,
                })
            } else {
                None
            }
        }

        fn decode(&self, input: u64) -> LitlenResult {
            let index = (input & self.primary_table_mask) as usize;
            let entry = self.primary_table[index];
            if entry & LITERAL_ENTRY != 0 {
                let input_bits = (entry & 0xf) as usize;
                let s1 = (entry >> 16) as u8;
                let s2 = (entry >> 24) as u8;

                let symbol_count = (entry & 0xf00) >> 8;
                match symbol_count {
                    1 => LitlenResult::SingleLiteral {
                        symbol: s1,
                        input_bits,
                    },
                    2 => LitlenResult::DoubleLiteral { s1, s2, input_bits },
                    _ => unreachable!(),
                }
            } else if entry & SECONDARY_TABLE_ENTRY != 0 {
                let input2 = input >> self.primary_table_bits;
                let index2 = (entry >> 16) + ((input2 as u32) & (entry & 0xff));
                let entry2 = self.secondary_table[index2 as usize];
                let input_bits = (entry2 & 0xf) as usize;
                let symbol = entry2 >> 4;
                LitlenResult::SecondaryTableLiteral { symbol, input_bits }
            } else {
                unreachable!("TODO: implement test covereage for this case")
            }
        }
    }

    #[test]
    fn test_rfc1951_example1() {
        // https://datatracker.ietf.org/doc/html/rfc1951 gives the following example
        // on page 8:
        //
        //    Symbol  Code
        //    ------  ----
        //    A       10
        //    B       0
        //    C       110
        //    D       111
        //
        // The code is completely defined by the sequence of bit lengths (2, 1, 3, 3).
        let t = LitlenTables::new(12, &[2, 1, 3, 3]).unwrap();
        assert_eq!(
            t.decode(0b_0_0_0000000_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 1,
                s2: 1,
                input_bits: 2
            },
        );
        assert_eq!(
            t.decode(0b_110_110_00_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 2,
                s2: 2,
                input_bits: 6
            },
        );
        assert_eq!(
            t.decode(0b_111_111_00_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 3,
                s2: 3,
                input_bits: 6
            },
        );
        assert_eq!(
            t.decode(0b_0_10_00000_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 1,
                s2: 0,
                input_bits: 3
            },
        );
    }

    #[test]
    fn test_rfc1951_example2() {
        // https://datatracker.ietf.org/doc/html/rfc1951 gives the following example
        // on page 9:
        //
        //    Symbol Length   Code
        //    ------ ------   ----
        //    A       3        010
        //    B       3        011
        //    C       3        100
        //    D       3        101
        //    E       3        110
        //    F       2         00
        //    G       4       1110
        //    H       4       1111
        let t = LitlenTables::new(12, &[3, 3, 3, 3, 3, 2, 4, 4]).unwrap();
        assert_eq!(
            t.decode(0b_010_011_00_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 0,
                s2: 1,
                input_bits: 6
            },
        );
        assert_eq!(
            t.decode(0b_00_00_0000_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 5,
                s2: 5,
                input_bits: 4
            },
        );
        assert_eq!(
            t.decode(0b_1111_1110_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 7,
                s2: 6,
                input_bits: 8
            },
        );
    }

    #[test]
    fn test_secondary_table() {
        // To smoke test the secondary table usage, we use a lopsided
        // tree that results in codes that are up to 15 bits long:
        //
        //    Symbol Length                 Code
        //    ------ ------   ------------------
        //    0       1                        0
        //    1       2                       10
        //    2       3                      110
        //    3       4                     1110
        //    4       5                   1_1110
        //    5       6                  11_1110
        //    6       7                 111_1110
        //    7       8                1111_1110
        //    8       9              1_1111_1110
        //    9       10            11_1111_1110
        //    10      11           111_1111_1110
        //    11      12          1111_1111_1110
        //    12      13        1_1111_1111_1110
        //    13      14       11_1111_1111_1110
        //    14      15      111_1111_1111_1110
        //    15      15      111_1111_1111_1111
        let t = LitlenTables::new(12, &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15])
            .unwrap();
        assert_eq!(
            t.decode(0b_0_0_000000_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 0,
                s2: 0,
                input_bits: 2
            },
        );
        assert_eq!(
            t.decode(0b_1110_1110_u8.reverse_bits() as u64),
            LitlenResult::DoubleLiteral {
                s1: 3,
                s2: 3,
                input_bits: 8
            },
        );
        assert_eq!(
            t.decode(0b_1111_1111_1111_1110u16.reverse_bits() as u64),
            LitlenResult::SecondaryTableLiteral {
                symbol: 15,
                input_bits: 15
            },
        );
        assert_eq!(
            t.decode(0b_1111_1111_1111_1111u16.reverse_bits() as u64),
            LitlenResult::SecondaryTableLiteral {
                symbol: 15,
                input_bits: 15
            },
        );
    }
}
