use std::convert::TryInto;

use simd_adler32::Adler32;

use crate::tables::{CLCL_ORDER, LEN_SYM_TO_LEN_BASE, LEN_SYM_TO_LEN_EXTRA};

const MAX_HEADER_BYTES: usize = 289;

#[derive(Debug)]
pub enum DecompressionError {
    /// Data cannot be decompressed with fdeflate, but may still be valid. The caller should
    /// fallback to a full zlib implementation.
    NotFDeflate,
    /// All input was consumed, but the end of the stream hasn't been reached.
    InsufficientInput,
    /// Too many literals were specified.
    InvalidHlit,
    /// The stream doesn't specify a valid huffman tree.
    BadHuffmanTree,
    /// The stream contains a distance code that was not allowed by the header.
    InvalidDistanceCode,
    /// The stream contains contains back-reference as the first symbol.
    InputStartsWithRun,
    /// The deflate stream checksum is incorrect.
    WrongChecksum,
}

#[derive(PartialEq, Eq)]
enum State {
    Header,
    Data,
    Checksum,
    Done,
}

/// Decompressor that reads fdeflate compressed streams.
pub struct Decompressor {
    buffer: u64,
    nbits: u8,
    bits_read: u64,
    data_table: [[u8; 2]; 4096],
    advance_table: [u16; 4096],

    state: State,

    // queued_input: Vec<u8>,
    queued_output: Option<(u8, usize)>,

    last: Option<u8>,

    checksum: Adler32,
}

impl Decompressor {
    /// Create a new decompressor.
    pub fn new() -> Self {
        Self {
            buffer: 0,
            nbits: 0,
            bits_read: 0,
            data_table: [[0; 2]; 4096],
            advance_table: [u16::MAX; 4096],
            queued_output: None,
            last: None,
            checksum: Adler32::new(),
            state: State::Header,
        }
    }

    fn fill_buffer(&mut self, input: &mut &[u8]) {
        if self.nbits == 64 {
            /* do nothing */
        } else if input.len() >= 8 {
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
        if input.len() < 5 {
            return Err(DecompressionError::InsufficientInput);
        }
        let mut block = &input[2..];

        // Validate zlib header follows PNG requirements.
        if input[0] & 0x0f != 0x08
            || (input[0] & 0xf0) > 0x70
            || u16::from_be_bytes(input[..2].try_into().unwrap()) % 31 != 0
        {
            return Err(DecompressionError::NotFDeflate);
        }

        // If FDICT is set, bail out and let the caller use a full zlib implementation.
        if input[1] & 0x20 != 0 {
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
            return Err(DecompressionError::InvalidHlit);
        }
        if hdist != 1 {
            return Err(DecompressionError::NotFDeflate);
        }

        // Read code length code lengths.
        let mut code_length_lengths = [0; 19];
        for i in 0..hclen {
            code_length_lengths[CLCL_ORDER[i]] =
                self.read_bits(3, &mut block)
                    .ok_or(DecompressionError::InsufficientInput)? as u8;
        }
        if code_length_lengths[16..] != [0; 3] || code_length_lengths[..16] != [4; 16] {
            return Err(DecompressionError::NotFDeflate);
        }
        let code_length_codes: [u16; 16] =
            crate::compute_codes(&code_length_lengths[..16].try_into().unwrap())
                .ok_or(DecompressionError::BadHuffmanTree)?;

        // Read literal/length code lengths.
        let mut lengths = [0; 286];
        for i in 0..hlit {
            let code = self
                .read_bits(4, &mut block)
                .ok_or(DecompressionError::InsufficientInput)? as u8;
            lengths[i] = code_length_codes[code as usize] as u8;
        }
        if lengths[0] == 0 || lengths.iter().any(|&l| l > 12) {
            return Err(DecompressionError::NotFDeflate);
        }

        // Read distance code lengths.
        let distance = code_length_codes[self
            .read_bits(4, &mut block)
            .ok_or(DecompressionError::InsufficientInput)?
            as usize];
        if distance != 1 {
            return Err(DecompressionError::NotFDeflate);
        }

        // Build the literal/length code table.
        let codes = crate::compute_codes(&lengths).ok_or(DecompressionError::BadHuffmanTree)?;

        // Check whether literal zero is assigned code zero. If so, our table can encode entries
        // with 3+ symbols even though each entry has only 2 data bytes.
        let use_extra_length = codes[0] == 0;

        for i in 0..256 {
            let code = codes[i];
            let length = lengths[i];
            let mut j = code;
            while j < 4096 && length != 0 {
                let extra_length = if use_extra_length {
                    ((j | 0xf000) >> length).trailing_zeros() as u8 / lengths[0]
                } else {
                    0
                };

                self.data_table[j as usize][0] = i as u8;
                self.advance_table[j as usize] =
                    (extra_length as u16 + 1) << 4 | (length + extra_length * lengths[0]) as u16;
                j += 1 << length;
            }

            if length > 0 && length <= 9 {
                for ii in 0..256 {
                    let code2 = codes[ii];
                    let length2 = lengths[ii];
                    if length2 != 0 && length + length2 <= 12 {
                        let mut j = code | (code2 << length);

                        while j < 4096 {
                            let extra_length = if use_extra_length {
                                ((j | 0xf000) >> (length + length2)).trailing_zeros() as u8
                                    / lengths[0]
                            } else {
                                0
                            };

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
                while j < 4096 && length != 0 {
                    self.advance_table[j as usize] = (i as u16 - 256) << 8 | (length as u16) << 4;
                    j += 1 << length;
                }
            }
        }

        let input_left = block.len();
        Ok(input.len() - input_left)
    }

    /// Decompresses a chunk of data.
    ///
    /// Returns the number of bytes read from `input` and the number of bytes written to `output`,
    /// or an error if the deflate stream is not valid. `input` is the compressed data. `output`
    /// is the buffer to write the decompressed data to. `end_of_input` indicates whether more
    /// data may be available in the future.
    ///
    /// If the stream is not in fdflate format, `NotFDeflate` is guaranteed to be returned *before
    /// consuming any input*. At this point, the same data should be passed to a full zlib
    /// implementation to decompress it.
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

        if let State::Done = self.state {
            return Ok((0, 0));
        }

        if let State::Header = self.state {
            if input.len() < MAX_HEADER_BYTES && !end_of_input {
                return Ok((0, 0));
            }
            let consumed = self.parse_header(input)?;
            remaining_input = &remaining_input[consumed..];
            self.state = State::Data;
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

        while let State::Data = self.state {
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
                    self.state = State::Checksum;
                    break;
                }

                let length_bits = LEN_SYM_TO_LEN_EXTRA[symbol - 257];
                let bits =
                    self.peak_bits(length_bits + advance_input_bits + 1) >> advance_input_bits;
                let length = LEN_SYM_TO_LEN_BASE[symbol - 257] + bits as usize;

                // Zero is the only valid distance code (corresponding to a distance of 1 byte).
                if (bits >> length_bits) != 0 {
                    return Err(DecompressionError::InvalidDistanceCode);
                }

                let last = if output_index > 0 {
                    output[output_index - 1]
                } else if let Some(last) = self.last {
                    last
                } else {
                    return Err(DecompressionError::InputStartsWithRun);
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

        if let State::Checksum = self.state {
            self.fill_buffer(&mut remaining_input);

            let align_bits = (8 - (self.bits_read % 8) as u8) % 8;
            if self.nbits >= 32 + align_bits {
                if align_bits != 0 {
                    self.consume_bits(align_bits);
                }
                #[cfg(not(fuzzing))]
                if (self.peak_bits(32) as u32).swap_bytes() != self.checksum.finish() {
                    return Err(DecompressionError::WrongChecksum);
                }
                self.state = State::Done;
                self.consume_bits(32);
            }
        }

        if self.state == State::Done || !end_of_input || output_index >= output.len() - 1 {
            if output_index > 0 {
                self.last = Some(output[output_index - 1])
            }
            let input_left = remaining_input.len();
            Ok((input.len() - input_left, output_index))
        } else {
            Err(DecompressionError::InsufficientInput)
        }
    }

    /// Returns true if the decompressor has finished decompressing the input.
    pub fn done(&self) -> bool {
        self.state == State::Done
    }
}

/// Decompresses the given input. Returns an error if the input is invalid or if the data is not
/// compressed with the fdeflate algorithm.
pub fn decompress_to_vec(input: &[u8]) -> Result<Vec<u8>, DecompressionError> {
    let mut decoder = Decompressor::new();
    let mut output = vec![0; 1024];
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
    use crate::tables::{LENGTH_TO_LEN_EXTRA, LENGTH_TO_SYMBOL};

    use super::*;
    use rand::Rng;

    fn roundtrip(data: &[u8]) {
        let compressed = crate::compress_to_vec(data);
        let decompressed = decompress_to_vec(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn tables() {
        for (i, &bits) in LEN_SYM_TO_LEN_EXTRA.iter().enumerate() {
            let len_base = LEN_SYM_TO_LEN_BASE[i];
            for j in 0..(1 << bits) {
                if i == 27 && j == 31 {
                    continue;
                }
                assert_eq!(LENGTH_TO_LEN_EXTRA[len_base + j - 3], bits, "{} {}", i, j);
                assert_eq!(
                    LENGTH_TO_SYMBOL[len_base + j - 3],
                    i as u16 + 257,
                    "{} {}",
                    i,
                    j
                );
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
