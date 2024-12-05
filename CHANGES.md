# Changelog

## 0.3.7

 - Optimized decoding of streams that use fixed huffman blocks ([#38], [#39])

[#38]: https://github.com/image-rs/fdeflate/pull/38
[#39]: https://github.com/image-rs/fdeflate/pull/39

## 0.3.6

This release improves end-to-end decoding peformance for PNG images by 10% on average,
with select images benefitting by as much as 50%.
These improvements were inspired by the algorithms in `zune-inflate`.

- Optimized building the Huffman table, which helps the performance of small images ([#31])
- Dropped the specialized fast path for a hardcoded Huffman table, which is no longer necessary ([#32])
- Add a fast path to the DEFLATE decoding loop that processes more bytes at a time, benefitting performance on large images ([#34])
- Re-add test files into the crates.io tarball since they are so small. They may be removed in the future if they grow in size ([#35])

[#31]: https://github.com/image-rs/fdeflate/pull/31
[#32]: https://github.com/image-rs/fdeflate/pull/32
[#34]: https://github.com/image-rs/fdeflate/pull/34
[#35]: https://github.com/image-rs/fdeflate/pull/35


## 0.3.5

- Fix handling of invalid inputs, so that errors are consistently detected
  regardless of how the input is "chunked" when feeding it into the
  decompressor.
- Add more fuzz testing.

## 0.3.4

- Fix bug where `Decompressor::read` might fail to return an error for a
  truncated deflate stream if the output buffer had room for exactly one more
  byte of data.

## 0.3.3

- Add `decompress_to_vec_bounded` method.

## 0.3.2

- Allow decoding into buffers without extra space.

## 0.3.1

- Strengthen postconditions on `Decompressor::read`.
- Add more fuzz testing.

## 0.3.0

- Added support for decoding arbitrary zlib streams.
