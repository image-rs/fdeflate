#[derive(Copy, Clone, Debug)]
pub enum Flush {
    Normal,
    Finish,
}

pub enum Status {
    Done,
    NeedMoreOutputSpace,
}

pub struct Compressor {}
impl Compressor {
    pub fn new() -> Compressor {
        Compressor {}
    }

    /// Returns (status, bytes_in, bytes_out)
    pub fn compress(
        &mut self,
        in_buf: &[u8],
        out_buf: &mut [u8],
        flush: Flush,
    ) -> (Status, usize, usize) {
        todo!()
    }
}

pub fn compress_to_vec(input: &[u8]) -> Vec<u8> {
    let mut compressor = Compressor::new();
    let mut output = vec![0; ::core::cmp::max(input.len() / 2, 2)];

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

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(data: &[u8]) {
        let compressed = compress_to_vec(data);
        let decompressed = miniz_oxide::inflate::decompress_to_vec(&compressed).unwrap();
        assert_eq!(&decompressed, data);
    }

    #[test]
    fn it_works() {
        roundtrip(b"Hello world!");
    }
}
