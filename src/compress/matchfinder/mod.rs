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

fn compute_hash3(v: u32) -> u32 {
    (0x330698ecu64.wrapping_mul(((v & 0xff_ffff) ^ 0x2722_0a95) as u64) >> 16) as u32
}

fn compute_hash(v: u64) -> u32 {
    // let mut hasher = fnv::FnvHasher::default();
    // std::hash::Hasher::write_u64(&mut hasher, v);
    // std::hash::Hasher::finish(&hasher) as u32

    (11400714785074694791u64.wrapping_mul(v) >> 40) as u32
}

fn compute_hash32(v: u32) -> u32 {
    // let mut hasher = fnv::FnvHasher::default();
    // std::hash::Hasher::write_u32(&mut hasher, v);
    // std::hash::Hasher::finish(&hasher) as u32

    (11400714785074694791u64.wrapping_mul(v as u64) >> 40) as u32
}

/// Find the length of the match between the current position and the previous position, searching
/// both forwards and backwards from the starting position.
fn match_length(
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
