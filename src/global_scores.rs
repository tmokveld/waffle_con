//! Incremental Myers bit vectors for unit-cost global edit distance to fixed reads.

use rustc_hash::FxHashSet;

const WORD_BITS: usize = usize::BITS as usize;
const HIGH_SHIFT: usize = WORD_BITS - 1;
const HIGH_BIT: usize = 1 << HIGH_SHIFT;
const ABSENT: u16 = u16::MAX;

struct ReadLayout {
    len: usize,
    words: usize,
    last_bit: usize,
}

/// Immutable equality masks shared by all branches of a search.
pub(crate) struct Patterns {
    reads: Vec<ReadLayout>,
    word_count: usize,
    symbol_rows: [u16; 256],
    // One contiguous row of words across all reads for each observed symbol.
    masks: Vec<usize>,
}

impl Patterns {
    /// `alphabet` contains every literal byte in `reads`.
    pub(crate) fn new(reads: &[&[u8]], alphabet: &FxHashSet<u8>) -> Self {
        let mut layouts = Vec::with_capacity(reads.len());
        let mut word_count = 0;
        for &read in reads {
            let len = read.len();
            let words = len.div_ceil(WORD_BITS);
            let last_bit = match len % WORD_BITS {
                // For an empty read this selects the initial boundary bit.
                0 => HIGH_BIT,
                bits => 1 << (bits - 1),
            };
            layouts.push(ReadLayout { len, words, last_bit });
            word_count += words;
        }

        let mut symbol_rows = [ABSENT; 256];
        for (row, &symbol) in (0u16..).zip(alphabet) {
            symbol_rows[usize::from(symbol)] = row;
        }
        let mut masks = vec![0; alphabet.len() * word_count];
        let mut read_start = 0;
        for (read, layout) in reads.iter().zip(&layouts) {
            for (i, &symbol) in read.iter().enumerate() {
                let row_start = usize::from(symbol_rows[usize::from(symbol)]) * word_count;
                masks[row_start + read_start + i / WORD_BITS] |= 1 << (i % WORD_BITS);
            }
            read_start += layout.words;
        }

        Self {
            reads: layouts,
            word_count,
            symbol_rows,
            masks,
        }
    }
}

/// The global DP columns for one candidate prefix, in input read order.
#[derive(Clone, Debug)]
pub(crate) struct Columns {
    read_count: usize,
    // Scores followed by interleaved positive and negative vertical delta words.
    state: Box<[usize]>,
}

impl Columns {
    pub(crate) fn root(patterns: &Patterns) -> Self {
        let read_count = patterns.reads.len();
        let mut state = Vec::with_capacity(read_count + 2 * patterns.word_count);
        state.extend(patterns.reads.iter().map(|read| read.len));
        for _ in 0..patterns.word_count {
            // D(i, 0) = i, so all vertical differences are +1.
            state.extend([usize::MAX, 0]);
        }
        Self {
            read_count,
            state: state.into_boxed_slice(),
        }
    }

    /// Append a literal byte. Use the same patterns as at the root.
    pub(crate) fn push(&mut self, patterns: &Patterns, symbol: u8) {
        let row = patterns.symbol_rows[usize::from(symbol)];
        if row == ABSENT {
            self.push_masks(patterns, std::iter::repeat(0));
        } else {
            let start = usize::from(row) * patterns.word_count;
            let masks = &patterns.masks[start..start + patterns.word_count];
            self.push_masks(patterns, masks.iter().copied());
        }
    }

    fn push_masks(&mut self, patterns: &Patterns, masks: impl Iterator<Item = usize>) {
        let (distances, vectors) = self.state.split_at_mut(self.read_count);
        let (vectors, _) = vectors.as_chunks_mut::<2>();
        let mut words = vectors.iter_mut().zip(masks);
        for (distance, read) in distances.iter_mut().zip(&patterns.reads) {
            // Seed each read with D(0, j) - D(0, j - 1) = +1 for global alignment.
            let mut ph = HIGH_BIT;
            let mut mh = 0;
            for (word, eq) in words.by_ref().take(read.words) {
                let ph_in = ph >> HIGH_SHIFT;
                let mh_in = mh >> HIGH_SHIFT;
                let pv = word[0];
                let mv = word[1];
                let xv = eq | mv;
                // Mh's high bit is also the carry out of (Eq & Pv) + Pv.
                // Setting Eq's low bit for Xh incorporates the incoming carry.
                let eq = eq | mh_in;
                let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
                ph = mv | !(xh | pv);
                mh = pv & xh;
                let shifted_ph = (ph << 1) | ph_in;
                let shifted_mh = (mh << 1) | mh_in;
                word[0] = shifted_mh | !(xv | shifted_ph);
                word[1] = xv & shifted_ph;
            }

            // Only the bit at the actual end of the read changes its score.
            *distance += (ph & read.last_bit != 0) as usize;
            *distance -= (mh & read.last_bit != 0) as usize;
        }
    }

    pub(crate) fn distances(&self) -> &[usize] {
        &self.state[..self.read_count]
    }
}

#[cfg(test)]
mod tests {
    use super::{Columns, Patterns};

    fn patterns_for(reads: &[&[u8]]) -> Patterns {
        let alphabet = reads.iter().flat_map(|&read| read).copied().collect();
        Patterns::new(reads, &alphabet)
    }

    /// The final row of scalar Wagner-Fischer gives one answer per prefix.
    fn scalar_prefix_distances(read: &[u8], candidate: &[u8]) -> Vec<usize> {
        let mut row: Vec<_> = (0..=candidate.len()).collect();
        for (i, &base) in read.iter().enumerate() {
            let mut diagonal = row[0];
            row[0] = i + 1;
            for (j, &symbol) in candidate.iter().enumerate() {
                let above = row[j + 1];
                row[j + 1] = (row[j] + 1)
                    .min(above + 1)
                    .min(diagonal + (base != symbol) as usize);
                diagonal = above;
            }
        }
        row
    }

    fn assert_prefixes(
        reads: &[&[u8]],
        patterns: &Patterns,
        columns: &mut Columns,
        candidate: &[u8],
        start: usize,
    ) {
        let expected: Vec<_> = reads
            .iter()
            .map(|read| scalar_prefix_distances(read, candidate))
            .collect();
        for len in start..=candidate.len() {
            if len > start {
                columns.push(patterns, candidate[len - 1]);
            }
            let distances = columns.distances();
            assert_eq!(distances.len(), reads.len());
            for (i, (&distance, row)) in distances.iter().zip(&expected).enumerate() {
                assert_eq!(distance, row[len], "read {i}, prefix length {len}");
            }
        }
    }

    fn assert_all_prefixes(reads: &[&[u8]], candidate: &[u8]) {
        let patterns = patterns_for(reads);
        let mut columns = Columns::root(&patterns);
        assert_prefixes(reads, &patterns, &mut columns, candidate, 0);
    }

    #[test]
    fn unequal_lengths_across_word_boundaries() {
        let alphabet = b"ACGTN*\0\xff";
        let template: Vec<_> = (0..151)
            .map(|i| alphabet[(i * 17 + i / 3 + i / 11) % alphabet.len()])
            .collect();
        let reads: &[&[u8]] = &[
            &template[7..7 + 65],
            b"",
            &template[4..4 + 128],
            &template[9..9 + 63],
            &template[..129],
            &template[70..71],
            &template[16..16 + 64],
            &template[8..8 + 127],
        ];
        let mut candidate = template.clone();
        candidate.insert(33, 0x80);
        candidate.remove(95);
        candidate[64] = b'C';
        assert_all_prefixes(reads, &candidate);
    }

    #[test]
    fn repeated_symbols_and_literal_ambiguities() {
        let repeated = [b'A'; 129];
        let reads: &[&[u8]] = &[&repeated, b"*N\0\xff", b"N*?NN"];
        let mut candidate = vec![b'A'; 132];
        candidate.extend_from_slice(b"\x80\x80N*?\0\xff");
        assert_all_prefixes(reads, &candidate);
    }

    #[test]
    fn all_byte_values() {
        let alphabet: Vec<_> = (0..=u8::MAX).collect();
        let reversed: Vec<_> = alphabet.iter().rev().copied().collect();
        assert_all_prefixes(&[&alphabet, b"\xff\0N*"], &reversed);
    }

    #[test]
    fn empty_baselines_and_no_reads() {
        assert_all_prefixes(&[b"", b""], b"A\0\xff");
        assert_all_prefixes(&[], b"A\0\xff");
    }

    #[test]
    fn cloned_prefixes_can_fork_independently() {
        let left_read = b"ACGT".repeat(33);
        let right_read = b"ACGTTGCA".repeat(17);
        let reads: &[&[u8]] = &[&left_read, b"", &right_read];
        let patterns = patterns_for(reads);
        let mut left = Columns::root(&patterns);
        let prefix = b"ACGT".repeat(19);
        assert_prefixes(reads, &patterns, &mut left, &prefix, 0);
        let mut right = left.clone();

        let mut left_candidate = prefix.clone();
        left_candidate.extend_from_slice(&left_read[prefix.len()..]);
        assert_prefixes(reads, &patterns, &mut left, &left_candidate, prefix.len());

        let mut right_candidate = prefix.clone();
        right_candidate.extend_from_slice(&right_read[prefix.len()..]);
        assert_prefixes(reads, &patterns, &mut right, &right_candidate, prefix.len());
    }
}
