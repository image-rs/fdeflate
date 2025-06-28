use super::compute_hash;

const CACHE_SIZE: usize = 1 << 18;

/// Find the length of the match between the current position and the previous position, searching
/// both forwards and backwards from the starting position.
fn match_length(data: &[u8], anchor: usize, mut ip: usize, mut prev_index: usize) -> (u16, usize) {
    assert!(
        prev_index < ip,
        "Match past current position: {prev_index} {ip}"
    );

    if data[ip] != data[prev_index] {
        return (0, ip);
    }

    let mut length = 1;
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

pub(crate) struct HashTableMatchFinder {
    hash_table: Box<[u32; CACHE_SIZE]>,

    /// Mask of low-bytes to consider for hashing.
    hash_mask: u64,
}
impl HashTableMatchFinder {
    pub(crate) fn new(min_match: u8) -> Self {
        assert!((3..=8).contains(&min_match));

        Self {
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
            hash_mask: if min_match == 8 {
                u64::MAX
            } else {
                (1 << (min_match.max(4) * 8)) - 1
            },
        }
    }

    pub(crate) fn get_and_insert(
        &mut self,
        data: &[u8],
        anchor: usize,
        ip: usize,
        value: u64,
        min_match: u16,
    ) -> (u16, u16, usize) {
        let min_offset = ip.saturating_sub(32768).max(1);

        let hash = compute_hash(value & self.hash_mask);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let offset = self.hash_table[hash_index] as usize;

        // Insert current value
        self.hash_table[hash_index] = ip as u32;

        if offset >= min_offset {
            let (length, start) = match_length(data, anchor, ip, offset);
            if length > min_match {
                return (length as u16, (ip - offset as usize) as u16, start);
            }
        }

        (0, 0, ip)
    }

    pub(crate) fn insert(&mut self, value: u64, offset: usize) {
        let hash = compute_hash(value & self.hash_mask);
        self.hash_table[(hash as usize) % CACHE_SIZE] = offset as u32;
    }
}
