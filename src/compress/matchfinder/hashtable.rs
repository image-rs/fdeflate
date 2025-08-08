use crate::compress::matchfinder::{Match, MatchFinder};

const CACHE_SIZE: usize = 1 << 16;

pub(crate) struct HashTableMatchFinder {
    hash_table: Box<[u32; CACHE_SIZE]>,
}
impl HashTableMatchFinder {
    pub(crate) fn new() -> Self {
        Self {
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
        }
    }
}
impl MatchFinder for HashTableMatchFinder {
    fn get_and_insert(
        &mut self,
        data: &[u8],
        base_index: u32,
        anchor: usize,
        ip: usize,
        value: u64,
    ) -> Match {
        let min_offset = (base_index + (ip as u32).saturating_sub(32768)).max(1);

        let hash = super::compute_hash(value);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let offset = self.hash_table[hash_index];

        // Insert current value
        self.hash_table[hash_index] = ip as u32 + base_index;

        if offset >= min_offset {
            assert!(
                ip > (offset - base_index) as usize,
                "ip={ip} offset={offset} base_index={base_index}"
            );
            let (length, start) = super::match_length::<true>(
                value,
                data,
                anchor,
                ip,
                (offset - base_index) as usize,
            );
            if length >= 8 {
                return Match::new(
                    length as u16,
                    (ip - (offset - base_index) as usize) as u16,
                    start,
                );
            }
        }

        Match::empty()
    }

    fn insert(&mut self, value: u64, index: u32) {
        let hash = super::compute_hash(value);
        self.hash_table[(hash as usize) % CACHE_SIZE] = index;
    }

    fn reset_indices(&mut self, old_base_index: u32) {
        for v in &mut *self.hash_table {
            *v = v.saturating_sub(old_base_index);
        }
    }
}
