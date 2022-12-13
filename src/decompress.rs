use std::convert::TryInto;

use simd_adler32::Adler32;

use crate::tables::{CLCL_ORDER, LEN_BITS, LEN_BASE};

const MAX_HEADER_BYTES: usize = 289;

#[derive(Debug)]
pub enum DecompressionError {
    /// Data cannot be decompressed with fdeflate, but may still be valid. The caller should
    /// fallback to a full zlib implementation.
    NotFDeflate,
    /// Data stream is corrupt.
    CorruptData,
}

pub struct Decompressor {
    buffer: u64,
    nbits: u8,
    bits_read: u64,
    data_table: [[u8; 2]; 4096],
    advance_table: [u16; 4096],

    read_header: bool,

    // queued_input: Vec<u8>,
    queued_output: Option<(u8, usize)>,

    done: bool,
    last: Option<u8>,

    checksum: Adler32,
}

impl Decompressor {
    pub fn new() -> Self {
        Self {
            buffer: 0,
            nbits: 0,
            bits_read: 0,
            data_table: [[0; 2]; 4096],
            advance_table: [u16::MAX; 4096],
            read_header: false,
            // queued_input: Vec::new(),
            queued_output: None,
            done: false,
            last: None,
            checksum: Adler32::new(),
        }
    }

    fn fill_buffer(&mut self, input: &mut &[u8]) {
        if input.len() >= 8 {
            self.buffer |= u64::from_le_bytes(input[..8].try_into().unwrap()) << self.nbits;
            *input = &mut &input[(63 - self.nbits as usize) / 8..];
            self.nbits |= 56;
        } else {
            let nbytes = input.len().min((64 - self.nbits as usize) / 8);
            let mut input_data = [0; 8];
            input_data[..nbytes].copy_from_slice(&input[..nbytes]);
            self.buffer |= u64::from_le_bytes(input_data) << self.nbits;
            self.nbits += nbytes as u8 * 8;
            *input = &mut &input[nbytes..];
        }
    }

    fn peak_bits(&mut self, nbits: u8) -> u64 {
        debug_assert!(nbits <= 56 && nbits <= self.nbits);
        self.buffer & ((1u64 << nbits) - 1)
    }
    fn consume_bits(&mut self, nbits: u8) {
        debug_assert!(self.nbits >= nbits);
        self.buffer >>= nbits;
        self.nbits -= nbits;
        self.bits_read += nbits as u64;
    }

    fn read_bits(&mut self, nbits: u8, input: &mut &[u8]) -> Option<u64> {
        if self.nbits < nbits {
            self.fill_buffer(input);
            if self.nbits < nbits {
                return None;
            }
        }

        let result = self.peak_bits(nbits);
        self.consume_bits(nbits);
        Some(result)
    }

    fn parse_header(&mut self, input: &[u8]) -> Result<usize, DecompressionError> {
        let mut block = &input[2..];

        // Validate zlib header follows PNG requirements.
        if input[0] & 0x0f != 0x08
            || (input[0] & 0xf0) > 0x70
            || u16::from_be_bytes(input[..2].try_into().unwrap()) % 31 != 0
        {
            return Err(DecompressionError::CorruptData);
        }

        // If FDICT is set, bail out and let the caller use a full zlib implementation.
        if block[1] & 0x10 != 0 {
            return Err(DecompressionError::NotFDeflate);
        }

        // fdeflate requires a single dynamic huffman block.
        let start = self.read_bits(3, &mut block).unwrap();
        if start != 0b101 {
            return Err(DecompressionError::NotFDeflate);
        }

        // Read the huffman table sizes.
        let hlit = self.read_bits(5, &mut block).unwrap() as usize + 257;
        let hdist = self.read_bits(5, &mut block).unwrap() as usize + 1;
        let hclen = self.read_bits(4, &mut block).unwrap() as usize + 4;
        if hlit > 286 {
            return Err(DecompressionError::CorruptData);
        }
        if hdist != 1 {
            return Err(DecompressionError::NotFDeflate);
        }

        // Read code length code lengths.
        let mut code_length_lengths = [0; 19];
        for i in 0..hclen {
            code_length_lengths[CLCL_ORDER[i]] = self.read_bits(3, &mut block).unwrap() as u8;
        }
        if code_length_lengths[16..] != [0; 3] || code_length_lengths[..16] != [4; 16] {
            return Err(DecompressionError::NotFDeflate);
        }
        let code_length_codes: [u16; 16] =
            crate::compute_codes(&code_length_lengths[..16].try_into().unwrap());

        // Read literal/length code lengths.
        let mut lengths = [0; 286];
        for i in 0..hlit {
            let code = self.read_bits(4, &mut block).unwrap() as u8;
            lengths[i] = code_length_codes[code as usize] as u8;
        }
        if lengths[0] == 0 || lengths.iter().any(|&l| l > 12) {
            return Err(DecompressionError::NotFDeflate);
        }

        // Read distance code lengths.
        let distance = code_length_codes[self.read_bits(4, &mut block).unwrap() as usize];
        if distance != 1 {
            return Err(DecompressionError::NotFDeflate);
        }

        // Build the literal/length code table.
        let codes = crate::compute_codes(&lengths);

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

        let input_left = block.len();
        Ok(input.len() - input_left)
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

        if self.done {
            return Ok((0, 0));
        }

        if !self.read_header {
            if input.len() < MAX_HEADER_BYTES && !end_of_input {
                return Ok((0, 0));
            }
            let consumed = self.parse_header(input)?;
            remaining_input = &remaining_input[consumed..];
            self.read_header = true;
        }

        if let Some((data, len)) = self.queued_output {
            let n = output.len().min(len);
            if data != 0 {
                output[..n].fill(data);
            }
            debug_assert_eq!(output_index, 0);
            output_index += n;
            self.queued_output = if n < len { Some((data, len - n)) } else { None };
        }

        loop {
            self.fill_buffer(&mut remaining_input);
            if self.nbits >= 48 {
                let bits = self.peak_bits(48);
                let advance0 = self.advance_table[(bits & 0xfff) as usize];
                let advance0_input_bits = advance0 & 0xf;
                let advance1 = self.advance_table[(bits >> advance0_input_bits) as usize & 0xfff];
                let advance1_input_bits = advance1 & 0xf;
                let advance2 = self.advance_table
                    [(bits >> (advance0_input_bits + advance1_input_bits)) as usize & 0xfff];
                let advance2_input_bits = advance2 & 0xf;
                let advance3 = self.advance_table[(bits
                    >> (advance0_input_bits + advance1_input_bits + advance2_input_bits))
                    as usize
                    & 0xfff];
                let advance3_input_bits = advance3 & 0xf;

                if advance0_input_bits > 0
                    && advance1_input_bits > 0
                    && advance2_input_bits > 0
                    && advance3_input_bits > 0
                {
                    let advance0_output_bytes = (advance0 >> 4) as usize;
                    let advance1_output_bytes = (advance1 >> 4) as usize;
                    let advance2_output_bytes = (advance2 >> 4) as usize;
                    let advance3_output_bytes = (advance3 >> 4) as usize;

                    if output_index
                        + advance0_output_bytes
                        + advance1_output_bytes
                        + advance2_output_bytes
                        + advance3_output_bytes
                        < output.len()
                    {
                        let data0 = self.data_table[(bits & 0xfff) as usize];
                        let data1 = self.data_table[(bits >> advance0_input_bits) as usize & 0xfff];
                        let data2 = self.data_table[(bits
                            >> (advance0_input_bits + advance1_input_bits))
                            as usize
                            & 0xfff];
                        let data3 = self.data_table[(bits
                            >> (advance0_input_bits + advance1_input_bits + advance2_input_bits))
                            as usize
                            & 0xfff];

                        let advance = advance0_input_bits
                            + advance1_input_bits
                            + advance2_input_bits
                            + advance3_input_bits;
                        self.consume_bits(advance as u8);

                        output[output_index] = data0[0];
                        output[output_index + 1] = data0[1];
                        output_index += advance0_output_bytes;
                        output[output_index] = data1[0];
                        output[output_index + 1] = data1[1];
                        output_index += advance1_output_bytes;
                        output[output_index] = data2[0];
                        output[output_index + 1] = data2[1];
                        output_index += advance2_output_bytes;
                        output[output_index] = data3[0];
                        output[output_index + 1] = data3[1];
                        output_index += advance3_output_bytes;
                        continue;
                    }
                }
            }

            if self.nbits < 18 {
                break;
            }

            let table_index = self.peak_bits(12);
            let data = self.data_table[table_index as usize];
            let advance = self.advance_table[table_index as usize];

            let advance_input_bits = (advance & 0x0f) as u8;
            let advance_output_bytes = (advance >> 4) as usize;

            if advance_input_bits > 0 {
                if output_index + 1 < output.len() {
                    output[output_index] = data[0];
                    output[output_index + 1] = data[1];
                    output_index += advance_output_bytes;
                    self.consume_bits(advance_input_bits);

                    if output_index > output.len() {
                        self.queued_output = Some((0, output_index - output.len()));
                        output_index = output.len();
                        break;
                    }
                } else if output_index + advance_output_bytes == output.len() {
                    debug_assert_eq!(advance_output_bytes, 1);
                    output[output_index] = data[0];
                    output_index += 1;
                    self.consume_bits(advance_input_bits);
                    break;
                } else {
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
                let bits =
                    self.peak_bits(length_bits + advance_input_bits + 1) >> advance_input_bits;
                let length = LEN_BASE[symbol - 257] + bits as usize;

                // Zero is the only valid distance code (corresponding to a distance of 1 byte).
                if (bits >> length_bits) != 0 {
                    return Err(DecompressionError::CorruptData);
                }

                let last = if output_index > 0 {
                    output[output_index - 1]
                } else if let Some(last) = self.last {
                    last
                } else {
                    return Err(DecompressionError::CorruptData);
                };

                self.consume_bits(length_bits + advance_input_bits + 1);

                // fdeflate only writes runs of zeros, but handling non-zero runs isn't hard and
                // it is too late to bail now.
                if last != 0 {
                    let end = (output_index + length).min(output.len());
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

        // self.data.extend_from_slice(&output[..output_index]);
        self.checksum.write(&output[..output_index]);
        // if self.done {
        //     if self.bits_read % 8 != 0 {
        //         self.consume_bits(8 - (self.bits_read % 8) as u8);
        //     }
        //     println!("current_nbits = {}", self.nbits);
        //     println!("remaining_input = {:x?}", remaining_input);
        //     if let Some(bits) = self.peak_bits(32, &mut remaining_input) {
        //         let full_checksum = simd_adler32::read::adler32(&mut std::io::Cursor::new(&*self.data)).unwrap();

        //         let checksum = self.checksum.finish();
        //         assert_eq!(checksum, full_checksum);
        //         assert_eq!(checksum.to_be(), bits as u32);
        //         if bits as u32 != checksum {
        //             panic!();
        //             return Err(DecompressionError::CorruptData);
        //         }
        //     }
        // }

        if self.done || !end_of_input || output_index >= output.len() - 1 {
            if output_index > 0 {
                self.last = Some(output[output_index - 1])
            }
            let input_left = remaining_input.len();
            Ok((input.len() - input_left, output_index))
        } else {
            Err(DecompressionError::CorruptData)
        }
    }

    pub fn done(&self) -> bool {
        self.done
    }
}

pub fn decompress_to_vec(input: &[u8]) -> Result<Vec<u8>, DecompressionError> {
    let mut decoder = Decompressor::new();
    let mut output = vec![0; 32 * 1024 * 1024];
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
    use crate::tables::{LEN_EXTRA, LEN_SYM};

    use super::*;
    use rand::Rng;

    fn roundtrip(data: &[u8]) {
        let compressed = crate::compress_to_vec(data);
        let decompressed = decompress_to_vec(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn tables() {
        for (i, &bits) in LEN_BITS.iter().enumerate() {
            let len_base = LEN_BASE[i];
            for j in 0..(1 << bits) {
                if i == 27 && j == 31 {
                    continue;
                }
                assert_eq!(LEN_EXTRA[len_base + j - 3], bits, "{} {}", i, j);
                assert_eq!(LEN_SYM[len_base + j - 3], i as u16 + 257, "{} {}", i, j);
            }
        }
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
