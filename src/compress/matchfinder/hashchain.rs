use crate::compress::{
    matchfinder::{Match, MatchFinder},
    WINDOW_SIZE,
};

const CACHE_SIZE: usize = 65536;

pub(crate) struct HashChainMatchFinder<const MIN_MATCH8: bool = false> {
    hash_table: Box<[u32; CACHE_SIZE]>,
    links: Box<[u32; WINDOW_SIZE]>,

    search_depth: u16,

    // /// If we already have a match of this length, limit lazy search to a smaller search depth.
    // good_length: u16,
    /// Stop searching for matches if the length is at least this long.
    nice_length: u16,

    min_match: u16,
    mask: u64,
}
impl<const MIN_MATCH8: bool> HashChainMatchFinder<MIN_MATCH8> {
    pub(crate) fn new(min_match: u16, search_depth: u16, nice_length: u16) -> Self {
        assert!(min_match >= 4 && min_match <= 8);
        assert_eq!(min_match == 8, MIN_MATCH8);

        Self {
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
            links: vec![0; WINDOW_SIZE].into_boxed_slice().try_into().unwrap(),
            search_depth,
            // good_length: 8,
            nice_length,
            min_match,
            mask: u64::MAX >> (8 * (8 - min_match)),
        }
    }
}

impl<const MIN_MATCH8: bool> MatchFinder for HashChainMatchFinder<MIN_MATCH8> {
    fn get_and_insert(
        &mut self,
        data: &[u8],
        base_index: u32,
        anchor: usize,
        ip: usize,
        value: u64,
    ) -> Match {
        let min_offset = (base_index + (ip as u32).saturating_sub(32768)).max(1);

        let mut best_offset = 0;
        let mut best_length = self.min_match - 1;
        let mut best_start = 0;

        let mut n = self.search_depth;
        // if min_match >= self.good_length {
        //     n >>= 2;
        // }

        let hash = super::compute_hash(value & self.mask);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let mut offset = self.hash_table[hash_index];

        // Insert current value
        self.hash_table[hash_index] = ip as u32 + base_index;
        self.links[ip % WINDOW_SIZE] = offset;

        // Visit previous matches
        loop {
            if offset < min_offset {
                break;
            }

            let (length, start) = super::match_length::<MIN_MATCH8>(
                value,
                data,
                anchor,
                ip,
                (offset - base_index) as usize,
            );
            if length > best_length {
                best_length = length;
                best_offset = offset;
                best_start = start;
            }
            if length >= self.nice_length || ip + length as usize == data.len() {
                break;
            }

            n -= 1;
            if n == 0 {
                break;
            }

            offset = self.links[offset as usize % WINDOW_SIZE];
        }

        if best_length >= self.min_match {
            return Match {
                length: best_length,
                distance: (ip - (best_offset - base_index) as usize) as u16,
                start: best_start,
            };
        }

        Match::empty()
    }

    fn insert(&mut self, value: u64, offset: u32) {
        let hash = super::compute_hash(value & self.mask);
        let prev_offset = self.hash_table[(hash as usize) % CACHE_SIZE];
        self.hash_table[(hash as usize) % CACHE_SIZE] = offset as u32;
        self.links[offset as usize % WINDOW_SIZE] = prev_offset;
    }

    fn reset_indices(&mut self, old_base_index: u32) {
        for v in &mut *self.hash_table {
            *v = v.saturating_sub(old_base_index);
        }

        for v in &mut *self.links {
            *v = v.saturating_sub(old_base_index);
        }
    }
}
