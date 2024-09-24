# Changelog

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
