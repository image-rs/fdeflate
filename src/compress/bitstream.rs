//! Methods for encoding the deflate bitstream.

use std::{
    collections::BinaryHeap,
    io::{self, Write},
};

use crate::{
    compress::BitWriter,
    tables::{
        BITMASKS, CLCL_ORDER, DIST_SYM_TO_DIST_BASE, DIST_SYM_TO_DIST_EXTRA, LENGTH_TO_LEN_EXTRA,
        LENGTH_TO_SYMBOL,
    },
};

pub(crate) fn distance_to_dist_sym(distance: u16) -> u8 {
    const LOOKUP: [u8; 16] = [0, 1, 2, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7];
    if distance <= 16 {
        return LOOKUP[distance as usize - 1];
    }

    let mut dist_sym = 29;
    while dist_sym > 0 && distance < DIST_SYM_TO_DIST_BASE[dist_sym as usize] {
        dist_sym -= 1;
    }
    dist_sym
}

pub(crate) enum Symbol {
    LiteralRun {
        start: u32,
        end: u32,
    },
    Backref {
        length: u16,
        distance: u16,
        dist_sym: u8,
    },
}

#[inline(always)]
pub(crate) fn write_block<W: Write>(
    writer: &mut BitWriter<W>,
    data: &[u8],
    base_index: u32,
    symbols: &[Symbol],
    eof: bool,
) -> io::Result<()> {
    let mut frequencies = [0u32; 286];
    let mut dist_frequencies = [0u32; 30];
    frequencies[256] = 1;

    let mut f2 = [0u32; 256];
    let mut f3 = [0u32; 256];
    let mut f4 = [0u32; 256];

    for symbol in symbols {
        match symbol {
            Symbol::LiteralRun { start, end } => {
                let mut chunks = data[(*start - base_index) as usize..(*end - base_index) as usize]
                    .chunks_exact(4);
                for chunk in &mut chunks {
                    frequencies[chunk[0] as usize] += 1;
                    f2[chunk[1] as usize] += 1;
                    f3[chunk[2] as usize] += 1;
                    f4[chunk[3] as usize] += 1;
                }
                for &lit in chunks.remainder() {
                    frequencies[lit as usize] += 1;
                }
            }
            Symbol::Backref {
                length, dist_sym, ..
            } => {
                let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                frequencies[sym] += 1;
                dist_frequencies[*dist_sym as usize] += 1;
            }
        }
    }

    for i in 0..256 {
        frequencies[i] += f2[i] + f3[i] + f4[i];
    }

    let mut lengths = [0u8; 286];
    let mut codes = [0u16; 286];
    build_huffman_tree(&frequencies, &mut lengths, &mut codes, 15);

    let mut dist_lengths = [0u8; 30];
    let mut dist_codes = [0u16; 30];
    build_huffman_tree(&dist_frequencies, &mut dist_lengths, &mut dist_codes, 15);

    let mut num_litlen_codes = 286;
    while num_litlen_codes > 257 && lengths[num_litlen_codes - 1] == 0 {
        num_litlen_codes -= 1;
    }

    let mut num_dist_codes = 30;
    while num_dist_codes > 1 && dist_lengths[num_dist_codes - 1] == 0 {
        num_dist_codes -= 1;
    }

    let mut code_length_frequencies = [0u32; 19];
    for &length in &lengths[..num_litlen_codes] {
        code_length_frequencies[length as usize] += 1;
    }
    for &length in &dist_lengths[..num_dist_codes] {
        code_length_frequencies[length as usize] += 1;
    }
    let mut code_length_lengths = [0u8; 19];
    let mut code_length_codes = [0u16; 19];
    build_huffman_tree(
        &code_length_frequencies,
        &mut code_length_lengths,
        &mut code_length_codes,
        7,
    );

    if eof {
        writer.write_bits(0b101, 3)?; // final block
    } else {
        writer.write_bits(0b100, 3)?; // non-final block
    }

    writer.write_bits(num_litlen_codes as u64 - 257, 5)?; // hlit
    writer.write_bits(num_dist_codes as u64 - 1, 5)?; // hdist
    writer.write_bits(15, 4)?; // hclen

    for j in 0..19 {
        writer.write_bits(code_length_lengths[CLCL_ORDER[j]] as u64, 3)?;
    }

    for &length in lengths[..num_litlen_codes]
        .iter()
        .chain(&dist_lengths[..num_dist_codes])
    {
        writer.write_bits(
            code_length_codes[length as usize] as u64,
            code_length_lengths[length as usize],
        )?;
    }

    for symbol in symbols {
        match symbol {
            Symbol::LiteralRun { start, end } => {
                let mut groups = data[(*start - base_index) as usize..(*end - base_index) as usize]
                    .chunks_exact(4);
                for group in &mut groups {
                    let code0 = codes[group[0] as usize] as u64;
                    let code1 = codes[group[1] as usize] as u64;
                    let code2 = codes[group[2] as usize] as u64;
                    let code3 = codes[group[3] as usize] as u64;

                    let len0 = lengths[group[0] as usize];
                    let len1 = lengths[group[1] as usize];
                    let len2 = lengths[group[2] as usize];
                    let len3 = lengths[group[3] as usize];

                    writer.write_bits(
                        code0
                            | (code1 << len0)
                            | (code2 << (len0 + len1))
                            | (code3 << (len0 + len1 + len2)),
                        len0 + len1 + len2 + len3,
                    )?;
                }

                for &lit in groups.remainder() {
                    writer.write_bits(codes[lit as usize] as u64, lengths[lit as usize])?;
                }
            }
            Symbol::Backref {
                length,
                distance,
                dist_sym,
            } => {
                let sym = LENGTH_TO_SYMBOL[*length as usize - 3] as usize;
                writer.write_bits(codes[sym] as u64, lengths[sym])?;
                let len_extra = LENGTH_TO_LEN_EXTRA[*length as usize - 3];
                let extra = (((*length as u32) - 3) & BITMASKS[len_extra as usize]) as u64;
                writer.write_bits(extra, len_extra)?;

                writer.write_bits(
                    dist_codes[*dist_sym as usize] as u64,
                    dist_lengths[*dist_sym as usize],
                )?;
                let dist_extra = DIST_SYM_TO_DIST_EXTRA[*dist_sym as usize];
                let extra = *distance - DIST_SYM_TO_DIST_BASE[*dist_sym as usize];

                writer.write_bits(extra as u64, dist_extra)?;
            }
        }
    }
    writer.write_bits(codes[256] as u64, lengths[256])?;
    Ok(())
}

fn build_huffman_tree(
    frequencies: &[u32],
    lengths: &mut [u8],
    codes: &mut [u16],
    length_limit: u8,
) -> bool {
    assert_eq!(frequencies.len(), lengths.len());
    assert_eq!(frequencies.len(), codes.len());

    if frequencies.iter().filter(|&&f| f > 0).count() <= 1 {
        lengths.fill(0);
        codes.fill(0);
        if let Some(i) = frequencies.iter().position(|&f| f > 0) {
            lengths[i] = 1;
        }
        return false;
    }

    #[derive(Eq, PartialEq, Copy, Clone, Debug)]
    struct Item(u32, u16);
    impl Ord for Item {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            other.0.cmp(&self.0)
        }
    }
    impl PartialOrd for Item {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    // Build a huffman tree
    let mut internal_nodes = Vec::new();
    let mut nodes = BinaryHeap::from_iter(
        frequencies
            .iter()
            .enumerate()
            .filter(|(_, &frequency)| frequency > 0)
            .map(|(i, &frequency)| Item(frequency, i as u16)),
    );
    while nodes.len() > 1 {
        let Item(frequency1, index1) = nodes.pop().unwrap();
        let mut root = nodes.peek_mut().unwrap();
        internal_nodes.push((index1, root.1));
        *root = Item(
            frequency1 + root.0,
            internal_nodes.len() as u16 + frequencies.len() as u16 - 1,
        );
    }

    // Walk the tree to assign code lengths
    lengths.fill(0);
    let mut stack = Vec::new();
    stack.push((nodes.pop().unwrap().1, 0));
    while let Some((node, depth)) = stack.pop() {
        let node = node as usize;
        if node < frequencies.len() {
            lengths[node] = depth as u8;
        } else {
            let (left, right) = internal_nodes[node - frequencies.len()];
            stack.push((left, depth + 1));
            stack.push((right, depth + 1));
        }
    }

    // Limit the codes to length length_limit
    let mut max_length = 0;
    for &length in lengths.iter() {
        max_length = max_length.max(length);
    }
    if max_length > length_limit {
        let mut counts = [0u32; 16];
        for &length in lengths.iter() {
            counts[length.min(length_limit) as usize] += 1;
        }

        let mut total = 0;
        for (i, count) in counts
            .iter()
            .enumerate()
            .skip(1)
            .take(length_limit as usize)
        {
            total += count << (length_limit as usize - i);
        }

        while total > 1u32 << length_limit {
            let mut i = length_limit as usize - 1;
            while counts[i] == 0 {
                i -= 1;
            }
            counts[i] -= 1;
            counts[length_limit as usize] -= 1;
            counts[i + 1] += 2;
            total -= 1;
        }

        // assign new lengths
        let mut len = length_limit;
        let mut indexes = frequencies.iter().copied().enumerate().collect::<Vec<_>>();
        indexes.sort_unstable_by_key(|&(_, frequency)| frequency);
        for &(i, frequency) in indexes.iter() {
            if frequency > 0 {
                while counts[len as usize] == 0 {
                    len -= 1;
                }
                lengths[i] = len;
                counts[len as usize] -= 1;
            }
        }
    }

    // Assign codes
    codes.fill(0);
    let mut code = 0u32;
    for len in 1..=length_limit {
        for (i, &length) in lengths.iter().enumerate() {
            if length == len {
                codes[i] = (code as u16).reverse_bits() >> (16 - len);
                code += 1;
            }
        }
        code <<= 1;
    }
    assert_eq!(code, 2 << length_limit);

    true
}
