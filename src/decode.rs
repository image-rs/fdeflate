const CLCL_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

const LEN_BITS: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];

const LEN_BASE: [usize; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 114, 131,
    163, 195, 227, 258,
];

const MAX_HEADER_BYTES: usize = 289;

#[derive(Debug)]
pub enum DecompressionError {
    /// Data cannot be decompressed with fdeflate, but may still be valid. The caller should
    /// fallback to a full zlib implementation.
    NotFDeflate(Vec<u8>),
    /// Data stream is corrupt.
    CorruptData,
}

pub struct Decompressor {
    buffer: u64,
    nbits: u8,
    data_table: [[u8; 2]; 4096],
    advance_table: [u16; 4096],

    header_contents: Option<Vec<u8>>,

    queued_input: Vec<u8>,
    queued_output: Option<(u8, usize)>,

    done: bool,
    last: Option<u8>,
}

impl Decompressor {
    pub fn new() -> Self {
        Self {
            buffer: 0,
            nbits: 0,
            data_table: [[0; 2]; 4096],
            advance_table: [u16::MAX; 4096],
            header_contents: Some(Vec::new()),
            queued_input: Vec::new(),
            queued_output: None,
            done: false,
            last: None,
        }
    }

    fn fill_buffer(&mut self, input: &[u8]) -> usize {
        let nbytes = input.len().min((64 - self.nbits as usize) / 8);

        let mut input_data = [0; 8];
        input_data[..nbytes].copy_from_slice(&input[..nbytes]);
        self.buffer |= u64::from_le_bytes(input_data) << self.nbits;
        self.nbits += nbytes as u8 * 8;
        nbytes
    }

    fn peak_bits(&mut self, nbits: u8, input: &mut &[u8]) -> Option<u64> {
        debug_assert!(nbits <= 56);

        if self.nbits < nbits {
            let advance = self.fill_buffer(input);
            *input = &input[advance..];
            if self.nbits < nbits {
                return None;
            }
        };

        Some(self.buffer & ((1u64 << nbits) - 1))
    }
    fn consume_bits(&mut self, nbits: u8) {
        debug_assert!(self.nbits >= nbits);
        self.buffer >>= nbits;
        self.nbits -= nbits;
    }

    fn read_bits(&mut self, nbits: u8, input: &mut &[u8]) -> Option<u64> {
        let result = self.peak_bits(nbits, input)?;
        self.consume_bits(nbits);
        Some(result)
    }

    fn parse_header(&mut self, header: Vec<u8>) -> Result<Vec<u8>, DecompressionError> {
        let mut block = &header[2..];

        // Validate zlib header follows PNG requirements.
        if header[0] & 0xf0 != 0x70
            || (header[0] & 0x0f) > 0x08
            || u16::from_be_bytes(header[..2].try_into().unwrap()) % 31 != 0
        {
            panic!();
            return Err(DecompressionError::CorruptData);
        }

        // If FDICT is set, bail out and let the caller use a full zlib implementation.
        if block[1] & 0x10 != 0 {
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }

        // fdeflate requires a single dynamic huffman block.
        let start = self.read_bits(3, &mut block).unwrap();
        if start != 0b101 {
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }

        // Read the huffman table sizes.
        let hlit = self.read_bits(5, &mut block).unwrap() as usize + 257;
        let hdist = self.read_bits(5, &mut block).unwrap() as usize + 1;
        let hclen = self.read_bits(4, &mut block).unwrap() as usize + 4;
        if hlit > 286 {
            return Err(DecompressionError::CorruptData);
        }
        if hdist != 1 {
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }

        // Read code length code lengths.
        let mut code_length_lengths = [0; 19];
        for i in 0..hclen {
            code_length_lengths[CLCL_ORDER[i]] = self.read_bits(3, &mut block).unwrap() as u8;
        }
        if code_length_lengths[16..] != [0; 3] || code_length_lengths[..16] != [4; 16] {
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }
        let code_length_codes: [u16; 16] =
            super::compute_codes(&code_length_lengths[..16].try_into().unwrap());

        // Read literal/length code lengths.
        let mut lengths = [0; 286];
        for i in 0..hlit {
            let code = self.read_bits(4, &mut block).unwrap() as u8;
            lengths[i] = code_length_codes[code as usize] as u8;
        }
        if lengths[0] == 0 || lengths.iter().any(|&l| l > 12){
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }

        // Read distance code lengths.
        let distance = code_length_codes[self.read_bits(4, &mut block).unwrap() as usize];
        if distance != 1 {
            panic!();
            return Err(DecompressionError::NotFDeflate(header));
        }

        // Build the literal/length code table.
        let codes = super::compute_codes(&lengths);

        for i in 0..256 {
            let code = codes[i];
            let length = lengths[i];
            let mut j = code;
            while j < 4096 {
                let extra_length = ((j | 0xf000) >> length).trailing_zeros() as u8 / lengths[0];

                self.data_table[j as usize][0] = i as u8;
                self.advance_table[j as usize] =
                    (extra_length as u16 + 1) << 4 | (length + extra_length * lengths[0]) as u16;
                j += 1 << length;
            }

            if length <= 9 {
                for ii in 0..256 {
                    let code2 = codes[ii];
                    let length2 = lengths[ii];
                    if length + length2 <= 12 {
                        let mut j = code | (code2 << length);
                        while j < 4096 {
                            let extra_length = ((j | 0xf000) >> (length + length2)).trailing_zeros()
                                as u8
                                / lengths[0];
                            self.data_table[j as usize][0] = i as u8;
                            self.data_table[j as usize][1] = ii as u8;
                            self.advance_table[j as usize] = (extra_length as u16 + 2) << 4
                                | (length + length2 + extra_length * lengths[0]) as u16;
                            j += 1 << (length + length2);
                        }
                    }
                }
            }
        }
        for i in 256..hlit {
            let code = codes[i];
            let length = lengths[i];
            if length != 0 {
                let mut j = code;
                while j < 4096 {
                    self.advance_table[j as usize] = (i as u16 - 256) << 8 | (length as u16) << 4;
                    j += 1 << length;
                }
            }
        }

        Ok(block.to_vec())
    }

    pub fn read(
        &mut self,
        input: &[u8],
        output: &mut [u8],
        end_of_input: bool,
    ) -> Result<(usize, usize), DecompressionError> {
        assert!(output.len() >= 2);
        debug_assert!(output.iter().all(|&b| b == 0));

        let mut remaining_input = &input[..];
        let mut output_index = 0;

        if let Some(mut header_contents) = self.header_contents.take() {
            let new_bytes = (MAX_HEADER_BYTES - header_contents.len()).min(remaining_input.len());
            header_contents.extend_from_slice(&remaining_input[..new_bytes]);
            remaining_input = &remaining_input[new_bytes..];

            if header_contents.len() < MAX_HEADER_BYTES && !end_of_input {
                // We need more bytes to parse the header.
                self.header_contents = Some(header_contents);
                return Ok((new_bytes, 0));
            }
            self.queued_input = self.parse_header(header_contents)?;
        }

        if let Some((data, len)) = self.queued_output {
            let n = output.len().min(len);
            if data != 0 {
                output[..n].fill(data);
            }
            debug_assert_eq!(output_index, 0);
            output_index += n;
            self.queued_output = if n < len { Some((data, len - n)) } else { None };
            if output_index == output.len() {
                return Ok((0, output_index));
            }
        }

        if !self.queued_input.is_empty() {
            let mut queued_input = std::mem::replace(&mut self.queued_input, Vec::new());
            let (consumed_in, consumed_out) =
                self.read(&queued_input, &mut output[output_index..], false)?;
            output_index += consumed_out;
            if consumed_in < queued_input.len() {
                queued_input.drain(..consumed_in);
                self.queued_input = queued_input;
                return Ok((0, output_index));
            }
        }

        while let Some(symbol) = self.peak_bits(12, &mut remaining_input) {
            let data = self.data_table[symbol as usize];
            let advance = self.advance_table[symbol as usize];

            let advance_input_bits = (advance & 0x0f) as u8;
            let advance_output_bytes = (advance >> 4) as usize;

            if advance_input_bits > 0 {
                if output_index >= output.len() {
                    break;
                } else if output_index + 1 == output.len() {
                    if advance_output_bytes == 1 {
                        output[output_index] = data[0];
                        output_index += advance_output_bytes;
                        self.consume_bits(advance_input_bits);
                    }
                    break;
                }

                output[output_index] = data[0];
                output[output_index + 1] = data[1];
                output_index += advance_output_bytes;
                self.consume_bits(advance_input_bits);

                if output_index > output.len() {
                    self.queued_output = Some((0, output_index - output.len()));
                    output_index = output.len();
                    break;
                }
            } else {
                let advance_input_bits = (advance_output_bytes & 0x0f) as u8;
                let symbol = 256 + (advance_output_bytes >> 4) as usize;

                // Check for end of input symbol.
                if symbol == 256 {
                    self.consume_bits(advance_input_bits);
                    self.done = true;
                    break;
                }

                let length_bits = LEN_BITS[symbol - 257];
                let bits = match self
                    .peak_bits(length_bits + advance_input_bits + 1, &mut remaining_input)
                {
                    Some(bits) => bits >> advance_input_bits,
                    None => break,
                };
                let length = LEN_BASE[symbol - 257] + bits as usize;

                // Zero is the only valid distance code (corresponding to a distance of 1 byte).
                if (bits >> length_bits) != 0 {
                    panic!();
                    return Err(DecompressionError::CorruptData);
                }

                let last = if output_index > 0 {
                    output[output_index - 1]
                } else if let Some(last) = self.last {
                    last
                } else {
                    panic!();
                    return Err(DecompressionError::CorruptData);
                };

                self.consume_bits(length_bits + advance_input_bits + 1);

                // fdeflate only writes runs of zeros, but handling non-zero runs isn't hard and
                // it is too late to bail now.
                if last != 0 {
                    let end = (output.len() + length).saturating_sub(output_index);
                    output[output_index..end].fill(last);
                }

                // The run can easily extend past the end of the output buffer. If so, queue the
                // output for the next call and break.
                if output_index + length > output.len() {
                    self.queued_output = Some((last, output_index + length - output.len()));
                    output_index = output.len();
                    break;
                } else {
                    output_index += length;
                }
            }
        }

        if self.done || !end_of_input || output_index >= output.len() - 1 {
            if output_index > 0 {
                self.last = Some(output[output_index - 1])
            }
            let input_left = remaining_input.len();
            Ok((input.len() - input_left, output_index))
        } else {
            panic!();
            Err(DecompressionError::CorruptData)
        }
    }

    pub fn done(&self) -> bool {
        self.done
    }
}

pub fn decompress_to_vec(input: &[u8]) -> Result<Vec<u8>, DecompressionError> {
    let mut decoder = Decompressor::new();
    let mut output = vec![0; 32 * 1024];
    let mut input_index = 0;
    let mut output_index = 0;
    while !decoder.done() {
        let (consumed, produced) =
            decoder.read(&input[input_index..], &mut output[output_index..], true)?;
        input_index += consumed;
        output_index += produced;
        output.resize(output_index + 32 * 1024, 0);
    }
    output.resize(output_index, 0);
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn roundtrip(data: &[u8]) {
        let compressed = crate::compress_to_vec(data);
        let decompressed = decompress_to_vec(&compressed).unwrap();
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

    #[bench]
    fn bench_runs(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = (rng.gen_range::<u8, _>(0..255)).saturating_sub(253);
        }
        let compressed = crate::compress_to_vec(&data);
        b.bytes = compressed.len() as u64;
        b.iter(|| decompress_to_vec(&compressed));
    }

    #[bench]
    fn bench_low(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = (rng.gen_range::<u8, _>(0..16)).wrapping_sub(8);
        }
        let compressed = crate::compress_to_vec(&data);
        b.bytes = compressed.len() as u64;
        b.iter(|| decompress_to_vec(&compressed));
    }

    #[bench]
    fn bench_mixture(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            if rng.gen_range(0..200) == 1 {
                *byte = rng.gen();
            } else {
                *byte = rng.gen_range::<u8, _>(0..32).wrapping_sub(16);
            }
        }
        let compressed = crate::compress_to_vec(&data);
        b.bytes = compressed.len() as u64;
        b.iter(|| decompress_to_vec(&compressed));
    }

    #[bench]
    fn bench_uniform_random(b: &mut test::Bencher) {
        let mut rng = rand::thread_rng();
        let mut data = vec![0; 1024 * 1024];
        for byte in &mut data {
            *byte = rng.gen();
        }
        let compressed = crate::compress_to_vec(&data);
        b.bytes = compressed.len() as u64;
        b.iter(|| decompress_to_vec(&compressed));
    }
}
