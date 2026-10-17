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

`StreamingCompressor` accepts input and output slices, returning separate consumed
and produced counts plus a stopping reason. It supports raw DEFLATE and zlib,
partial/sync/full flushes, and finishing. Pass unconsumed input again with the
same flush request. Once all input has been accepted, that flush remains pending
until delivered; repeating a completed flush without new input is a no-op.
`process_uninit` writes to `MaybeUninit<u8>` storage without initializing unused
capacity, using only safe Rust.

Working storage is bounded independently of the stream length. The encoder
buffers at most 128 KiB of new input plus a 32 KiB history window, limits symbol
storage to 16,388 entries, and reserves 131,201 bytes for encoded output. Each
block compares fixed, dynamic, and stored encoding costs; the stored alternative
bounds expansion. A sufficiently large destination receives a batch directly;
smaller destinations drain the reserved buffer before more input is accepted.
The match tables and Huffman workspace have fixed bounds too. Automatic blocks
preserve history, and long streams rebase match indices without resetting it.
`reset` reuses the working allocations. These block boundaries can change the
compressed bytes and compression ratio compared with the previous writer API.

`Compressor<W>` drives the same core through a reusable 32 KiB writer buffer.
After an underlying writer error it must be discarded. For resumable output
backpressure, use `StreamingCompressor` and retain any bytes not yet written.

`Decompressor::read_progress` reports consumed and produced counts even on an
error. Its output prefix must contain the previous 32 KiB of history (or all
previous output if shorter), as with `read`. A decoding error remains terminal
until `reset`.

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

The algorithms in this crate take inspiration many other compression libraries including:
[fpnge](https://github.com/veluca93/fpnge),
[lz4](https://github.com/lz4/lz4),
[xz-utils](https://tukaani.org/xz/)
[zlib-ng](https://github.com/zlib-ng/zlib-ng),
[zlib-rs](https://github.com/trifectatechfoundation/zlib-rs),
[ZStandard](https://github.com/facebook/zstd),
and
[zune-inflate](https://github.com/etemesi254/zune-image/tree/main/zune-inflate).
