mod hashchain;
mod hashtable;

pub(crate) use hashchain::HashChainMatchFinder;
pub(crate) use hashtable::HashTableMatchFinder;

pub(crate) struct Match {
    pub length: u16,
    pub distance: u16,
    pub start: usize,
}
impl Match {
    pub fn new(length: u16, distance: u16, match_start: usize) -> Self {
        debug_assert!(length >= 3, "Match length must be at least 3");
        Self {
            length,
            distance,
            start: match_start,
        }
    }

    pub fn empty() -> Self {
        Self {
            length: 0,
            distance: 0,
            start: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn end(&self) -> usize {
        self.start + self.length as usize
    }
}

pub(crate) trait MatchFinder {
    fn get_and_insert(
        &mut self,
        data: &[u8],
        base_index: u32,
        anchor: usize,
        ip: usize,
        value: u64,
    ) -> Match;

    fn insert(&mut self, value: u64, offset: u32);

    fn reset_indices(&mut self, old_base_index: u32);
}

fn compute_hash(v: u64) -> u32 {
    (11400714785074694791u64.wrapping_mul(v) >> 40) as u32
}

/// Find the length of the match between the current position and the previous position.
///
/// This function checks that `value` matches the 4 bytes (or 8 bytes if MIN_MATCH >= 8) at the
/// previous index, and if so, searches both forwards and backwards to find the starting position
/// and total length of the match.
fn match_length<const MIN_MATCH8: bool>(
    value: u64,
    data: &[u8],
    anchor: usize,
    mut ip: usize,
    mut prev_index: usize,
) -> (u16, usize) {
    assert!(
        prev_index < ip,
        "Match past current position: {prev_index} {ip}"
    );

    let prev = u64::from_ne_bytes(data[prev_index..][..8].try_into().unwrap());

    let mut length;
    if MIN_MATCH8 {
        if value != prev {
            return (0, ip);
        }
        length = 8;
    } else {
        if value as u32 != prev as u32 {
            return (0, ip);
        }
        length = (value ^ prev).trailing_zeros() as usize / 8;
    }

    // Search backwards to find where the match starts.
    while length < 258 && ip > anchor && prev_index > 0 && data[ip - 1] == data[prev_index - 1] {
        length += 1;
        ip -= 1;
        prev_index -= 1;
    }

    // Search forwards to find the full length of the match.
    'fsearch: {
        let slice_length = (data.len() - ip - length).min(258 - length);
        let mut chunks = data[ip + length..][..slice_length].chunks_exact(8);
        let mut prev_chunks = data[prev_index + length..][..slice_length].chunks_exact(8);

        for (chunk, prev_chunk) in (&mut chunks).zip(&mut prev_chunks) {
            if chunk == prev_chunk {
                length += 8;
            } else {
                let chunk = u64::from_ne_bytes(chunk.try_into().unwrap());
                let prev_chunk = u64::from_ne_bytes(prev_chunk.try_into().unwrap());
                length += (chunk ^ prev_chunk).trailing_zeros() as usize / 8;
                break 'fsearch; // skip the remainder loop below
            }
        }
        for (chunk, prev_chunk) in chunks.remainder().iter().zip(prev_chunks.remainder()) {
            if *chunk != *prev_chunk {
                break;
            }
            length += 1;
        }
    }

    (length as u16, ip)
}
