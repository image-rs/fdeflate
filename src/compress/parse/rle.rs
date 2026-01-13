use std::io::{self, Write};

use crate::compress::{matchfinder::NullMatchFinder, parse::ParserInner, BitWriter, Flush};

pub(crate) struct RleParser {
    inner: ParserInner<NullMatchFinder>,
}

impl RleParser {
    pub fn new(skip_ahead_shift: u8) -> Self {
        Self {
            inner: ParserInner::new(skip_ahead_shift, NullMatchFinder),
        }
    }

    pub fn reset_indices(&mut self, old_base_index: u32) {
        self.inner.reset_indices(old_base_index);
    }

    pub fn compress<W: Write>(
        &mut self,
        writer: &mut BitWriter<W>,
        data: &[u8],
        base_index: u32,
        start: usize,
        flush: Flush,
    ) -> io::Result<usize> {
        self.inner.start_compress(data, base_index, start);

        let lookahead = if flush == Flush::None { 258 } else { 7 };
        let max_ip = data.len().saturating_sub(lookahead);

        loop {
            let m = self.inner.advance_to_match(data, base_index, max_ip);
            if m.is_empty() {
                break;
            }

            self.inner.ip = m.end();
            self.inner.insert_match(base_index, m);
            self.inner
                .write_block_if_ready(writer, data, base_index, flush)?;
        }

        self.inner
            .end_compress(writer, data, base_index, start, flush)
    }
}
