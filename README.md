# fdeflate

[![crates.io](https://img.shields.io/crates/v/fdeflate.svg)](https://crates.io/crates/fdeflate)
[![Documentation](https://docs.rs/fdeflate/badge.svg)](https://docs.rs/fdeflate)
[![Build Status](https://img.shields.io/github/actions/workflow/status/image-rs/fdeflate/rust.yml?label=Rust%20CI)](https://github.com/image-rs/fdeflate/actions)

A fast, safe, and modern zlib implementation for PNG.

## Overview

This crate contains an optimized implementation of the [deflate algorithm](https://en.wikipedia.org/wiki/Deflate) tuned for PNG images, but capable of handling arbitrary data. It was created to serve as the compression backend for the [`png`](https://crates.io/crates/png) crate and can also be used as a standalone library. The upcoming 0.4.x series will feature a fully compatible streaming encoder supporting uncompressed output, traditional levels 1-9, and two specialized compression levels designed for PNG images.

## Decompression

The decoder employs several modern features than help it achieve exceptional performance that rivals or exceeds the performance of `zlib-ng` and `zlib-rs` without using any `unsafe` code:

* It uses larger Huffman tables with up to 4096 entries, to better make use of the larger L1 caches available on modern processors.
* [Multi-byte literal decoding](https://fastcompression.blogspot.com/2015/10/huffman-revisited-part-4-multi-bytes.html) enables decoding two literals with a single table lookup. Since PNG files often have a higher percentage of literals, this is particularly helpful for them.
* Optimized Huffman table construction which requires fewer random writes.

## Compression

The compressor is still under active development, with the current status tracked on the main branch.

<!-- Internally, `fdeflate` (like other compression libraries) maps different compression levels to different algorithms and parameters. -->

[Preliminary benchmark data](https://github.com/image-rs/fdeflate/discussions/68) shows that `fdeflate` meaningfully outperforms `zlib-rs` for levels 1-3 and is slightly better for levels 4-7. Levels 8 and 9 aren't yet implemented.

### Ultra-fast compression

The "ultra-fast" compression level is a specialized compression mode that is several _times_ faster than other modes while still providing some amount of compression. It is designed to be performance-competitive with [`QOI`](https://qoiformat.org/) while still being compatible with zlib. It does so by making a bunch of simplifying assumptions:

- Exactly one block per deflate stream.
- No distance codes except for run length encoding of zeros.
- A single fixed Huffman tree trained on a large corpus of PNG images.
- All Huffman codes are <= 12 bits.

In the 0.3.x series, the ultra-fast mode is the only supported compression mode.

## Inspiration

The algorithms in this crate take inspiration many other compression including libraries:
[fpnge](https://github.com/veluca93/fpnge),
[lz4](https://github.com/lz4/lz4),
[xz-utils](https://tukaani.org/xz/)
[zlib-ng](https://github.com/zlib-ng/zlib-ng),
[zlib-rs](https://github.com/trifectatechfoundation/zlib-rs),
[ZStandard](https://github.com/facebook/zstd),
and
[zune-inflate](https://github.com/etemesi254/zune-image/tree/main/zune-inflate).
