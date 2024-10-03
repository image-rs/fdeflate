use crate::{
    decompress::{EXCEPTIONAL_ENTRY, LITERAL_ENTRY, SECONDARY_TABLE_ENTRY},
    tables::{LEN_SYM_TO_LEN_BASE, LEN_SYM_TO_LEN_EXTRA},
    DecompressionError,
};

fn next_codeword(mut codeword: u16, table_size: u16) -> u16 {
    if codeword == table_size - 1 {
        return codeword;
    }

    // If the codeword isn't all ones (which is the final code), then we need to
    // advance it to the next code.
    //
    // This bit manipulation is equivelent to `(codeword.reverse_bits() + 1).reverse_bits()`
    // but is more efficient.
    let adv = (u16::BITS - 1) - (codeword ^ (table_size - 1)).leading_zeros();
    let bit = 1 << adv;
    codeword &= bit - 1;
    codeword |= bit;
    codeword
}

pub fn build_table(
    lengths: &[u8],
    entries: &[u32],
    codes: &mut [u16],
    primary_table: &mut [u32],
    secondary_table: &mut Vec<u16>,
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
    let mut sorted_symbols = [0; 288];
    for symbol in 0..lengths.len() {
        let length = lengths[symbol];
        sorted_symbols[offsets[length as usize]] = symbol;
        offsets[length as usize] += 1;
    }

    let mut codeword = 0u16;
    let mut i = histogram[0];

    let primary_table_bits = primary_table.len().ilog2() as usize;

    // Populate the primary decoding table
    for length in 1..=primary_table_bits {
        let current_table_end = 1 << length;

        // Loop over all symbols with the current code length and set their table entries.
        for _ in 0..histogram[length] {
            let symbol = sorted_symbols[i];
            i += 1;

            primary_table[codeword as usize] =
                entries.get(symbol).cloned().unwrap_or((symbol as u32) << 4) | length as u32;

            codes[symbol] = codeword;
            codeword = next_codeword(codeword, current_table_end as u16);
        }

        // If we aren't at the maximum table size, double the size of the table.
        if length < primary_table_bits {
            primary_table.copy_within(0..current_table_end, current_table_end);
        }
    }

    let secondary_table_bits = max_length - primary_table_bits;
    let secondary_table_size = 1 << secondary_table_bits;
    let secondary_table_mask = secondary_table_size - 1;

    let mut subtable_start = 0;
    let mut subtable_prefix = !0;
    secondary_table.clear();
    for length in (primary_table_bits + 1)..=max_length {
        for _ in 0..histogram[length] {
            if codeword & 0xfff != subtable_prefix {
                subtable_prefix = codeword & 0xfff;
                subtable_start = secondary_table.len();
                primary_table[subtable_prefix as usize] = ((subtable_start as u32) << 16)
                    | EXCEPTIONAL_ENTRY
                    | SECONDARY_TABLE_ENTRY
                    | secondary_table_mask as u32;
                secondary_table.resize(subtable_start + secondary_table_size, 0);
            }

            let symbol = sorted_symbols[i];
            i += 1;

            codes[symbol] = codeword;

            let mut s = codeword >> primary_table_bits;
            while s < secondary_table_size as u16 {
                debug_assert_eq!(secondary_table[subtable_start + s as usize], 0);
                secondary_table[subtable_start + s as usize] =
                    ((symbol as u16) << 4) | (length as u16);
                s += 1 << (length - primary_table_bits);
            }

            // compression.secondary_table[subtable_start + (codeword >> 12) as usize] =
            //     ((symbol as u16) << 4) | length as u16;

            codeword = next_codeword(codeword, 1 << length);
        }

        // if length < max_length && codeword & 0xfff == subtable_prefix {
        //     compression.litlen_table[subtable_prefix] = subtable_start as u32
        //         | ((length - 12 + 1) << 8) as u32
        //         | EXCEPTIONAL_ENTRY
        //         | SECONDARY_TABLE_ENTRY;

        //     compression
        //         .secondary_table
        //         .extend_from_within(subtable_start..);
        // }
    }

    true
}
