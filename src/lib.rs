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

struct Table {
    code: [u16; 286],
    nbits: [u8; 286],
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

        // Write code lengths for literal/length alphabet
        for _ in 0..255 {
            self.write_bits(0b0001, 4, &mut out_buf); // 8 bits for first 255 codes
        }
        for _ in 0..2 {
            self.write_bits(0b1001, 4, &mut out_buf); // 9 bits for last 2 codes
        }

        // Write code lengths for distance alphabet
        for _ in 0..1 {
            self.write_bits(0, 4, &mut out_buf); // 5 bits per distance code
        }

        // Write data
        for &byte in in_buf {
            if byte != 255 {
                self.write_bits(byte.reverse_bits() as u64, 8, &mut out_buf);
            } else {
                self.write_bits(0b1_1111_1110, 9, &mut out_buf);
            }
        }

        // Write end of block
        self.write_bits(0b1_1111_1111, 9, &mut out_buf);
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

    fn roundtrip(data: &[u8]) {
        let compressed = compress_to_vec(data);
        let decompressed = miniz_oxide::inflate::decompress_to_vec_zlib(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn it_works() {
        roundtrip(b"Hello world!");
    }
}
