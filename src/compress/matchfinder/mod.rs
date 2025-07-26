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

/// Find the 4+ byte length of the match between the current position and the previous position.
///
/// This function checks that `value` matches the 4 bytes at the previous index, and if so, searches
/// both forwards and backwards to find the starting position and total length of the match.
fn match_length4(
    value: u32,
    data: &[u8],
    anchor: usize,
    mut ip: usize,
    mut prev_index: usize,
) -> (u16, usize) {
    assert!(
        prev_index < ip,
        "Match past current position: {prev_index} {ip}"
    );

    if value != u32::from_ne_bytes(data[prev_index..][..4].try_into().unwrap()) {
        return (0, ip);
    }

    let mut length = 4;
    while length < 258 && ip > anchor && prev_index > 0 && data[ip - 1] == data[prev_index - 1] {
        length += 1;
        ip -= 1;
        prev_index -= 1;
    }
    while length < 258 && ip + length < data.len() && data[ip + length] == data[prev_index + length]
    {
        length += 1;
    }
    (length as u16, ip)
}

/// Find the 8+ byte length of the match between the current position and the previous position.
///
/// This function checks that `value` matches the 8 bytes at the previous index, and if so, searches
/// both forwards and backwards to find the starting position and total length of the match.
fn match_length8(
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

    if value != u64::from_ne_bytes(data[prev_index..][..8].try_into().unwrap()) {
        return (0, ip);
    }

    let mut length = 8;
    while length < 258 && ip > anchor && prev_index > 0 && data[ip - 1] == data[prev_index - 1] {
        length += 1;
        ip -= 1;
        prev_index -= 1;
    }
    while length < 258 && ip + length < data.len() && data[ip + length] == data[prev_index + length]
    {
        length += 1;
    }
    (length as u16, ip)
}
