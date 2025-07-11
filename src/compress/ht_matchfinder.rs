use super::compute_hash;

const CACHE_SIZE: usize = 1 << 16;

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

pub(crate) struct HashTableMatchFinder {
    hash_table: Box<[u32; CACHE_SIZE]>,
}
impl HashTableMatchFinder {
    pub(crate) fn new() -> Self {
        Self {
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
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

        let hash = compute_hash(value);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let offset = self.hash_table[hash_index] as usize;

        // Insert current value
        self.hash_table[hash_index] = ip as u32;

        if offset >= min_offset {
            let (length, start) = match_length(value, data, anchor, ip, offset);
            if length > min_match {
                return (length as u16, (ip - offset as usize) as u16, start);
            }
        }

        (0, 0, ip)
    }

    pub(crate) fn insert(&mut self, value: u64, offset: usize) {
        let hash = compute_hash(value);
        self.hash_table[(hash as usize) % CACHE_SIZE] = offset as u32;
    }
}
