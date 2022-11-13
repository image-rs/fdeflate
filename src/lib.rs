use std::io::Write;

#[derive(Copy, Clone, Debug)]
pub enum Flush {
    Normal,
    Finish,
}

pub enum Status {
    Done,
    NeedMoreOutputSpace,
}

// Dynamic programming huffman tree algorithm from fpnge.
fn compute_code_lengths(
    freqs: &[u64],
    min_limit: &[u8],
    max_limit: &[u8],
    calculated_nbits: &mut [u8],
) {
    debug_assert_eq!(freqs.len(), min_limit.len());
    debug_assert_eq!(freqs.len(), max_limit.len());
    debug_assert_eq!(freqs.len(), calculated_nbits.len());
    let len = freqs.len();

    for i in 0..len {
        debug_assert!(min_limit[i] >= 1);
        debug_assert!(min_limit[i] <= max_limit[i]);
    }

    let precision = *max_limit.iter().max().unwrap();
    let num_patterns = 1 << precision;

    let mut dynp = vec![u64::MAX; (num_patterns + 1) * (len + 1)];
    let index = |sym: usize, off: usize| sym * (num_patterns + 1) + off;

    dynp[index(0, 0)] = 0;
    for sym in 0..len {
        for bits in min_limit[sym]..=max_limit[sym] {
            let off_delta = 1 << (precision - bits);
            for off in 0..=num_patterns.saturating_sub(off_delta) {
                dynp[index(sym + 1, off + off_delta)] = dynp[index(sym, off)]
                    .saturating_add(freqs[sym] * u64::from(bits))
                    .min(dynp[index(sym + 1, off + off_delta)]);
            }
        }
    }

    let mut sym = len;
    let mut off = num_patterns;

    while sym > 0 {
        sym -= 1;
        assert!(off > 0);

        for bits in min_limit[sym]..=max_limit[sym] {
            let off_delta = 1 << (precision - bits);
            if off_delta <= off
                && dynp[index(sym + 1, off)]
                    == dynp[index(sym, off - off_delta)].saturating_add(freqs[sym] * u64::from(bits))
            {
                off -= off_delta;
                calculated_nbits[sym] = bits;
                break;
            }
        }
    }

    for i in 0..len {
        debug_assert!(calculated_nbits[i] >= min_limit[i]);
        debug_assert!(calculated_nbits[i] <= max_limit[i]);
    }
}

fn build_literal_length_tree(buf: &[u8]) -> Vec<u8> {
    let mut counts = vec![1u64; 47];
    counts[46] = 1;
    for &b in buf {
        match b {
            0..=15 => counts[b as usize] += 1,
            240..=255 => counts[b as usize - 240 + 16] += 1,
            b => counts[((b - 16) >> 4) as usize + 32] += 1,
        }
    }

    let mut min_limit = vec![1u8; 47];
    let max_limit = vec![8u8; 47];
    for i in 32..46 {
        min_limit[i] = 8;
    }

    let mut lengths = vec![0u8; 47];
    compute_code_lengths(&counts, &min_limit, &max_limit, &mut lengths);

    let mut output_lengths = Vec::with_capacity(257);
    output_lengths.extend_from_slice(&lengths[..16]);
    for &len in &lengths[32..46] {
        output_lengths.extend_from_slice(&[len + 4; 16]);
    }
    output_lengths.extend_from_slice(&lengths[16..32]);
    output_lengths.push(lengths[46]);
    assert_eq!(output_lengths.len(), 257);

    for i in 16..240 {
        assert_eq!(output_lengths[i], 12);
    }

    output_lengths
}

fn compute_codes(lengths: &[u8]) -> Vec<u16> {
    let mut codes = vec![0u16; lengths.len()];

    let max_length = *lengths.iter().max().unwrap();

    let mut code = 0u16;
    for len in 1..=max_length {
        for (i, &length) in lengths.iter().enumerate() {
            if length == len {
                codes[i] = code.reverse_bits() >> (16 - len);
                code += 1;
            }
        }
        code <<= 1;
    }

    assert_eq!(code, 2 << max_length);

    codes
}

pub struct Compressor {
    buffer: u64,
    nbits: u8,
    bytes_written: usize,
}
impl Compressor {
    pub fn new() -> Compressor {
        Compressor {
            buffer: 0,
            nbits: 0,
            bytes_written: 0,
        }
    }

    fn write_bits(&mut self, bits: u64, nbits: u8, output: &mut &mut [u8]) {
        debug_assert!(nbits <= 64);

        self.buffer |= bits << self.nbits;
        self.nbits += nbits;

        if self.nbits >= 64 {
            output.write_all(&self.buffer.to_le_bytes()).unwrap();
            self.nbits -= 64;
            self.buffer = bits >> (nbits - self.nbits);
            self.bytes_written += 8;
        }
    }

    fn flush_to_byte_boundary(&mut self, output: &mut &mut [u8]) {
        if self.nbits % 8 != 0 {
            self.write_bits(0, 8 - self.nbits % 8, output);
        }
    }

    fn flush(&mut self, output: &mut &mut [u8]) {
        self.flush_to_byte_boundary(output);
        if self.nbits > 0 {
            output
                .write_all(&self.buffer.to_le_bytes()[..self.nbits as usize / 8])
                .unwrap();
            self.bytes_written += self.nbits as usize / 8;
            self.buffer = 0;
            self.nbits = 0;
        }
    }

    /// Returns (status, bytes_in, bytes_out)
    pub fn compress(
        &mut self,
        in_buf: &[u8],
        mut out_buf: &mut [u8],
        flush: Flush,
    ) -> (Status, usize, usize) {
        self.write_bits(0x0178, 16, &mut out_buf); // zlib header

        self.write_bits(0b1, 1, &mut out_buf); // BFINAL
        self.write_bits(0b10, 2, &mut out_buf); // Dynamic Huffman block

        self.write_bits(0, 5, &mut out_buf); // 257 length / literal codes
        self.write_bits(0, 5, &mut out_buf); // 1 distance code
        self.write_bits(15, 4, &mut out_buf); // 19 code length codes

        // Write code lengths for code length alphabet
        for _ in 0..3 {
            self.write_bits(0, 3, &mut out_buf);
        }
        for _ in 0..16 {
            self.write_bits(4, 3, &mut out_buf);
        }

        let lengths = build_literal_length_tree(in_buf);
        assert_eq!(lengths.len(), 257);
        for &len in &lengths {
            self.write_bits(
                (len.reverse_bits() >> 4) as u64,
                4,
                &mut out_buf,
            );
        }

        // Write code lengths for distance alphabet
        for _ in 0..1 {
            self.write_bits(0, 4, &mut out_buf); // 5 bits per distance code
        }

        // Write data
        let codes = compute_codes(&lengths);
        for &byte in in_buf {
            self.write_bits(
                codes[byte as usize] as u64,
                lengths[byte as usize],
                &mut out_buf,
            );
        }

        // Write end of block
        self.write_bits(codes[256] as u64, lengths[256], &mut out_buf);
        self.flush(&mut out_buf);

        // Write Adler32 checksum
        let checksum: u32 = adler32::adler32(in_buf).unwrap();
        out_buf.write_all(checksum.to_be_bytes().as_ref()).unwrap();
        self.bytes_written += 4;

        (Status::Done, in_buf.len(), self.bytes_written)
    }
}

pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new();
    let mut output = vec![0; ::core::cmp::max(1024 + input.len() * 2, 2)];

    let mut in_pos = 0;
    let mut out_pos = 0;
    loop {
        let (status, bytes_in, bytes_out) =
            compressor.compress(&input[in_pos..], &mut output[out_pos..], Flush::Finish);

        out_pos += bytes_out;
        in_pos += bytes_in;

        match status {
            Status::Done => {
                output.truncate(out_pos);
                break;
            }
            Status::NeedMoreOutputSpace => {
                // We need more space, so resize the vector.
                if output.len().saturating_sub(out_pos) < 30 {
                    output.resize(output.len() * 2, 0)
                }
            }
        }
    }

    output.resize(out_pos, 0);
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn roundtrip(data: &[u8]) {
        let compressed = compress_to_vec(data);
        let decompressed = miniz_oxide::inflate::decompress_to_vec_zlib(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn it_works() {
        roundtrip(b"Hello world!");
    }

    #[test]
    fn constant() {
        roundtrip(&vec![0; 2048]);
        roundtrip(&vec![5; 2048]);
        roundtrip(&vec![128; 2048]);
        roundtrip(&vec![254; 2048]);
    }

    #[test]
    fn random() {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 2048];
        for _ in 0..10 {
            for byte in &mut data {
                *byte = rng.gen();
            }
            roundtrip(&data);
        }
    }
}
