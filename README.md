# fdeflate

[![crates.io](https://img.shields.io/crates/v/image.svg)](https://crates.io/crates/fdeflate)
[![Documentation](https://docs.rs/image/badge.svg)](https://docs.rs/fdeflate)
[![Build Status](https://github.com/image-rs/image/workflows/Rust%20CI/badge.svg)](https://github.com/image-rs/fdeflate/actions)

A fast deflate implementation.

This crate contains an optimized implementation of the deflate algorithm tuned to compress PNG
images. It is compatible with standard zlib, but make a bunch of simplifying assumptions that
drastically improve encoding performance:

- Exactly one block per deflate stream.
- No distance codes except for run length encoding of zeros.
- A single fixed huffman tree trained on a large corpus of PNG images.
- All huffman codes are <= 12 bits\*.

\*This restriction may be removed in the future.

