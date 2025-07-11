use super::{compute_hash, compute_hash3, WINDOW_SIZE};

const CACHE3_SIZE: usize = 1 << 15;
const CACHE_SIZE: usize = 1 << 16;

/// Find the length of the match between the current position and the previous position, searching
/// both forwards and backwards from the starting position.
fn match_length(data: &[u8], ip: usize, prev_index: usize) -> u16 {
    assert!(
        prev_index < ip,
        "Match past current position: {prev_index} {ip}"
    );

    let mut length = 0;
    while length < 258 && ip + length < data.len() && data[ip + length] == data[prev_index + length]
    {
        length += 1;
    }
    length as u16
}

fn left_child(index: usize) -> usize {
    2 * (index as usize % WINDOW_SIZE)
}

fn right_child(index: usize) -> usize {
    2 * (index as usize % WINDOW_SIZE) + 1
}

/// Match finder that uses a binary tree to find matches.
///
/// Based on bt_matchfinder.h from libdeflate.
pub(crate) struct BTreeMatchFinder {
    hash3_table: Option<Box<[u32; CACHE3_SIZE]>>,
    hash_table: Box<[u32; CACHE_SIZE]>,
    child_links: Box<[u32; WINDOW_SIZE * 2]>,
    search_depth: u16,
    early_return_length: usize,
}
impl BTreeMatchFinder {
    pub(crate) fn new(min_match: u8) -> Self {
        assert!((3..=4).contains(&min_match));

        Self {
            hash3_table: (min_match == 3)
                .then(|| vec![0; CACHE3_SIZE].into_boxed_slice().try_into().unwrap()),
            hash_table: vec![0; CACHE_SIZE].into_boxed_slice().try_into().unwrap(),
            child_links: vec![0; WINDOW_SIZE * 2]
                .into_boxed_slice()
                .try_into()
                .unwrap(),
            search_depth: 2000,
            early_return_length: 256,
        }
    }

    fn update(
        &mut self,
        data: &[u8],
        ip: usize,
        value: u64,
        min_match: u16,
        record_matches: bool,
    ) -> (u16, u16, usize) {
        let min_offset = ip.saturating_sub(WINDOW_SIZE).max(1);

        let mut best_offset = 0;
        let mut best_length = min_match - 1;

        // Handle 3-byte matches
        if let Some(hash3_table) = &mut self.hash3_table {
            let hash3 = compute_hash3(value as u32);
            if best_length < min_match && min_match <= 3 {
                let hash3_offset = hash3_table[(hash3 as usize) % CACHE3_SIZE] as usize;
                if hash3_offset >= ip.saturating_sub(8192).max(1) {
                    let length = match_length(data, ip, hash3_offset);
                    if length >= 3 {
                        best_length = length;
                        best_offset = hash3_offset as u32;
                    }
                }
            }
            hash3_table[(hash3 as usize) % CACHE3_SIZE] = ip as u32;
        }

        // Lookup current value
        let hash = compute_hash(value & 0xffff_ffff);
        let hash_index = (hash as usize) % CACHE_SIZE;
        let mut offset = self.hash_table[hash_index] as usize;
        self.hash_table[hash_index] = ip as u32;

        let mut pending_left = left_child(ip);
        let mut pending_right = right_child(ip);

        if offset < min_offset {
            self.child_links[pending_left] = 0;
            self.child_links[pending_right] = 0;
            return (0, 0, ip);
        }

        let mut best_left_length = 0;
        let mut best_right_length = 0;
        let mut length = 0;

        // Visit previous matches
        // eprintln!("---");
        let mut depth_remaining = self.search_depth;
        loop {
            if data[ip + length] == data[offset + length] {
                while length < 258
                    && ip + length < data.len()
                    && data[ip + length] == data[offset + length]
                {
                    length += 1;
                }

                // for i in 0..length.min(self.early_return_length) {
                //     assert_eq!(
                //         data[ip + i],
                //         data[offset + i],
                //         "{i} {length} ip={ip} data_len={}",
                //         data.len()
                //     );
                // }

                if record_matches && length > best_length as usize {
                    best_length = length as u16;
                    best_offset = offset as u32;
                }

                if length >= self.early_return_length || ip + length == data.len() {
                    self.child_links[pending_left] = self.child_links[left_child(offset)];
                    self.child_links[pending_right] = self.child_links[right_child(offset)];
                    break;
                }
            }

            assert!(ip + length < data.len());

            if data[offset + length] < data[ip + length] {
                self.child_links[pending_left] = offset as u32;
                pending_left = right_child(offset);
                offset = self.child_links[pending_left] as usize;

                best_left_length = length;
                if best_right_length < length {
                    length = best_right_length;
                }
                // length = length.min(best_right_length);
                // eprintln!(
                //     "left {best_right_length},{best_left_length} dist={}",
                //     ip - offset
                // );
            } else {
                assert!(
                    data[offset + length] > data[ip + length],
                    "{length} {depth_remaining} {offset} {min_offset}"
                );

                self.child_links[pending_right] = offset as u32;
                pending_right = left_child(offset);
                offset = self.child_links[pending_right] as usize;

                best_right_length = length;
                if best_left_length < length {
                    length = best_left_length;
                }
                // length = length.min(best_left_length);
                // eprintln!(
                //     "right {best_right_length},{best_left_length} dist={}",
                //     ip - offset
                // );
            }

            depth_remaining -= 1;
            if offset <= min_offset || depth_remaining == 0 {
                self.child_links[pending_left] = 0;
                self.child_links[pending_right] = 0;
                break;
            }
        }

        if best_length >= min_match {
            return (best_length as u16, (ip - best_offset as usize) as u16, ip);
        }

        (0, 0, ip)
    }

    pub(crate) fn get_and_insert(
        &mut self,
        data: &[u8],
        ip: usize,
        value: u64,
        min_match: u16,
    ) -> (u16, u16, usize) {
        self.update(data, ip, value, min_match, true)
    }

    pub(crate) fn insert(&mut self, data: &[u8], value: u64, ip: usize) {
        self.update(data, ip, value, 3, false);
    }
}
