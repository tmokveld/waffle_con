use rustc_hash::FxHashMap;

// Counts are reset through touched entries; maps retain fresh-map growth order.
// Ordering depends on std HashMap internals and is guarded against real fresh maps.
pub(crate) struct CandidateScratch {
    counts: [usize; 256],
    touched: [u8; 256],
    touched_len: usize,
    levels: [FxHashMap<u8, usize>; 8],
    level: usize,
}

impl Default for CandidateScratch {
    fn default() -> Self {
        Self {
            counts: [0; 256],
            touched: [0; 256],
            touched_len: 0,
            levels: std::array::from_fn(|_| FxHashMap::default()),
            level: 0,
        }
    }
}

impl CandidateScratch {
    pub(crate) fn reset(&mut self) {
        for &symbol in &self.touched[..self.touched_len] {
            self.counts[symbol as usize] = 0;
        }
        if self.touched_len >= 2 {
            self.levels[self.level].clear();
        }
        self.touched_len = 0;
        self.level = 0;
    }

    pub(crate) fn record(&mut self, symbol: u8) {
        let index = symbol as usize;
        if self.counts[index] != 0 {
            self.counts[index] += 1;
            return;
        }

        self.touched[self.touched_len] = symbol;
        self.touched_len += 1;
        self.counts[index] = 1;
        if self.touched_len == 1 {
            return;
        }
        if self.levels[self.level].capacity() == 0 {
            self.levels[self.level].reserve(1);
        }
        if self.touched_len == 2 {
            self.levels[0].entry(self.touched[0]).or_insert(0);
        }
        if self.levels[self.level].len() == self.levels[self.level].capacity() {
            let next_level = self.level + 1;
            let (earlier, later) = self.levels.split_at_mut(next_level);
            let previous = &mut earlier[self.level];
            let next = &mut later[0];
            if next.capacity() == 0 {
                next.reserve(previous.capacity() + 1);
            }
            // Fresh-map growth replays occupied buckets in this same traversal order.
            for &key in previous.keys() {
                next.entry(key).or_insert(0);
            }
            previous.clear();
            self.level = next_level;
        }
        self.levels[self.level].entry(symbol).or_insert(0);
    }

    pub(crate) fn ordered_counts(&self) -> impl Iterator<Item = (u8, usize)> + '_ {
        let singleton = (self.touched_len == 1)
            .then(|| (self.touched[0], self.counts[self.touched[0] as usize]));
        singleton.into_iter().chain(
            self.levels[self.level]
                .keys()
                .map(|&symbol| (symbol, self.counts[symbol as usize])),
        )
    }

    /// First-seen membership, not fresh-map traversal order.
    pub(crate) fn distinct_symbols(&self) -> &[u8] {
        &self.touched[..self.touched_len]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compare_stream(scratch: &mut CandidateScratch, stream: &[u8]) {
        scratch.reset();
        let mut fresh = FxHashMap::<u8, usize>::default();
        assert!(scratch
            .ordered_counts()
            .eq(fresh.iter().map(|(&b, &n)| (b, n))));
        for &symbol in stream {
            *fresh.entry(symbol).or_insert(0) += 1;
            scratch.record(symbol);
            assert_eq!(
                scratch.ordered_counts().collect::<Vec<_>>(),
                fresh.iter().map(|(&b, &n)| (b, n)).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn test_fresh_map_order_across_growth_and_duplicates() {
        let mut scratch = CandidateScratch::default();
        for reverse in [false, true] {
            let mut bytes: Vec<u8> = (0..=255).collect();
            if reverse {
                bytes.reverse();
            }
            let mut stream = Vec::new();
            for &symbol in &bytes {
                if let Some(&previous) = stream.last() {
                    stream.push(previous);
                }
                stream.extend_from_slice(&[symbol, symbol]);
            }
            compare_stream(&mut scratch, &stream);
        }
    }

    #[test]
    fn test_large_small_empty_reset_order() {
        let mut scratch = CandidateScratch::default();
        let full: Vec<u8> = (0..=255)
            .flat_map(|b| std::iter::repeat(b).take(1 + b as usize % 5))
            .collect();
        let streams = [
            vec![],
            vec![255, 255],
            full,
            vec![0, 255, 0],
            vec![],
            vec![0, 4, 8, 255, 4],
        ];
        for _ in 0..3 {
            for stream in &streams {
                compare_stream(&mut scratch, stream);
            }
            for stream in streams.iter().rev() {
                compare_stream(&mut scratch, stream);
            }
        }
    }

    #[test]
    fn test_singleton_ambiguity_transitions_match_fresh_maps() {
        let mut scratch = CandidateScratch::default();
        let full: Vec<u8> = (0..=255)
            .flat_map(|byte| std::iter::repeat(byte).take(1 + byte as usize % 5))
            .collect();
        let streams = [
            vec![],
            vec![0],
            vec![0, 0, 0],
            vec![0, 0, 255, 0, 255],
            full,
            vec![255, 255],
            vec![],
            vec![4, 4, 8, 4],
            vec![4],
            vec![],
        ];
        for _ in 0..3 {
            for stream in &streams {
                compare_stream(&mut scratch, stream);
            }
            for stream in streams.iter().rev() {
                compare_stream(&mut scratch, stream);
            }
        }
    }
}
