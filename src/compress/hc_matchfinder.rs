use super::{compute_hash, compute_hash3, WINDOW_SIZE};

const CACHE3_SIZE: usize = 1 << 15;
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

pub(crate) struct HashChainMatchFinder {
    hash3_table: Option<Box<[u32; CACHE3_SIZE]>>,
    hash_table: Box<[u32; CACHE_SIZE]>,
    links: Box<[u32; WINDOW_SIZE]>,

    search_depth: u16,

    // /// If we already have a match of this length, limit lazy search to a smaller search depth.
    // good_length: u16,

    /// Stop searching for matches if the length is at least this long.
    nice_length: u16,

    /// Mask of low-bytes to consider for hashing.
    hash_mask: u64,
}
impl HashChainMatchFinder {
    pub(crate) fn new(search_depth: u16, nice_length: u16, min_match: u8) -> Self {
        assert!((3..=8).contains(&min_match));

        Self {
            hash3_table: (min_match == 3)
                .then(|| vec![0; CACHE3_SIZE].into_boxed_slice().try_into().unwrap()),
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
            links: vec![0; WINDOW_SIZE].into_boxed_slice().try_into().unwrap(),
            search_depth,
            // good_length: 8,
            nice_length,
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

        let mut best_offset = 0;
        let mut best_length = min_match - 1;
        let mut best_ip = 0;

        let mut n = self.search_depth;
        // if min_match >= self.good_length {
        //     n >>= 2;
        // }

        let hash = compute_hash(value & self.hash_mask);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let mut offset = self.hash_table[hash_index] as usize;

        // Insert current value
        self.hash_table[hash_index] = ip as u32;
        self.links[ip % WINDOW_SIZE] = offset as u32;

        // Visit previous matches
        loop {
            if offset < min_offset {
                break;
            }

            let (length, start) = match_length(data, anchor, ip, offset);
            if length > best_length {
                best_length = length;
                best_offset = offset as u32;
                best_ip = start;
                // } else if best_length > min_match {
                //     break;
            }
            if length >= self.nice_length || ip + length as usize == data.len() {
                break;
            }

            n -= 1;
            if n == 0 {
                break;
            }

            offset = self.links[offset % WINDOW_SIZE] as usize;
        }

        if let Some(hash3_table) = &mut self.hash3_table {
            let hash3 = compute_hash3(value as u32);
            if best_length < min_match && min_match <= 3 {
                let hash3_offset = hash3_table[(hash3 as usize) % CACHE3_SIZE] as usize;
                if hash3_offset >= ip.saturating_sub(8192).max(1) {
                    let (length, start) = match_length(data, anchor, ip, hash3_offset);
                    if length >= 3 {
                        best_length = length;
                        best_offset = hash3_offset as u32;
                        best_ip = start;
                    }
                }
            }
            hash3_table[(hash3 as usize) % CACHE3_SIZE] = ip as u32;
        }

        if best_length >= min_match {
            return (
                best_length as u16,
                (ip - best_offset as usize) as u16,
                best_ip,
            );
        }

        (0, 0, ip)
    }

    pub(crate) fn insert(&mut self, value: u64, offset: usize) {
        if let Some(hash3_table) = &mut self.hash3_table {
            let hash3 = compute_hash3(value as u32);
            hash3_table[(hash3 as usize) % CACHE3_SIZE] = offset as u32;
        }

        let hash = compute_hash(value & self.hash_mask);
        let prev_offset = self.hash_table[(hash as usize) % CACHE_SIZE];
        self.hash_table[(hash as usize) % CACHE_SIZE] = offset as u32;
        self.links[offset as usize % WINDOW_SIZE] = prev_offset;
    }
}
