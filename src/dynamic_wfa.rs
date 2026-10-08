
use rustc_hash::FxHashMap as HashMap;
use simple_error::bail;
use crate::candidate_scratch::CandidateScratch;

/// Lifecycle of a [`DWFALite`].
/// A DWFA starts [`Active`](DWFALiteState::Active). 
/// Finalizing it or passing the edit-distance cap moves it out of that state.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub enum DWFALiteState {
    /// The DWFA can still be extended.
    #[default]
    Active,
    /// `finalize` has completed and the DWFA can no longer be extended.
    Finalized,
    /// Another edit would pass `max_edit_distance`. The distance stays at the cap and this DWFA stops voting.
    ExceededEditDistanceLimit,
}

/// Fixed options for a [`DWFALite`].
/// These values are chosen when the DWFA is constructed and do not change as it is extended.
///
/// ```
/// use waffle_con::dynamic_wfa::{DWFALite, DWFALiteConfigBuilder};
/// let config = DWFALiteConfigBuilder::default()
///     .wildcard(Some(b'N'))
///     .allow_early_termination(true)
///     .max_edit_distance(Some(500))
///     .build()
///     .unwrap();
/// let dwfa = DWFALite::new(config);
/// assert_eq!(dwfa.max_edit_distance(), Some(500));
/// ```
#[derive(derive_builder::Builder, Clone, Debug, Eq, Hash, PartialEq)]
#[builder(default)]
pub struct DWFALiteConfig {
    /// Optional wildcard symbol that matches anything.
    pub wildcard: Option<u8>,
    /// When true, reaching the end of the baseline does not penalize a longer other sequence.
    pub allow_early_termination: bool,
    /// Absolute edit-distance cap. `None` means unlimited.
    pub max_edit_distance: Option<usize>,
}

// explicitly set our default for now, even if it's derivable
#[allow(clippy::derivable_impls)]
impl Default for DWFALiteConfig {
    fn default() -> Self {
        DWFALiteConfig {
            wildcard: None,
            allow_early_termination: false,
            max_edit_distance: None,
        }
    }
}

/// The core dynamic WFA structure for a lite implementation.
/// It is structured such that all the sequences being built are maintained **outside** of this struct (hence the "lite").
/// Essentially, it is keep the wavefront information, but all sequences are updated and tracked elsewhere.
/// As a result, all updates must pass the sequences so the information can get tracked.
/// If the outside sequences are changed in any way other than appending, then this may become desynchronized.
/// Conceptually, if this is a 2D grid, the baseline sequence will go from top to bottom (y-axis) and the other sequence will go from left to right (x-axis).
/// This means each character we add to `other_seq` will add a new _column_ to the grid.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct DWFALite {
    /// Fixed matching and edit-distance options for this DWFA.
    config: DWFALiteConfig,
    /// Whether this DWFA is still active, finalized, or frozen at its edit-distance cap.
    state: DWFALiteState,
    /// The minimum edit distance from `other_seq` to `baseline_seq` allowing for a large port of the tail of the `baseline_seq` to get ignored.
    edit_distance: usize,
    /// This is the core wavefront for the current `edit_distance`.
    /// It is always length 2*`edit_distance`+1.
    /// For our mental model, the median diagonal means we have used the same number of bases in our baseline and other.
    /// Anything below the median diagonal means we have used more baseline than other.
    /// Anything above the median diagonal means we have used more other than baseline.
    /// The value stored will always be the number of bases consumed in `other_seq`.
    /// After an active update without early termination, at least one equals the length of `other_seq` - `offset`.
    wavefront: Vec<usize>,
    // this is an offset into `other_seq` that the baseline starts
    offset: usize,
}

impl Default for DWFALite {
    fn default() -> Self {
        DWFALite {
            config: DWFALiteConfig::default(),
            state: DWFALiteState::Active,
            edit_distance: 0,
            wavefront: vec![0],
            offset: 0,
        }
    }
}
    

impl DWFALite {
    /// Creates a new dynamic WFA instance from fixed options.
    /// # Arguments
    /// * `config` - wildcard, early termination, and the absolute edit-distance cap
    pub fn new(config: DWFALiteConfig) -> DWFALite {
        DWFALite {
            config,
            ..Default::default()
        }
    }

    /// Sets an offset into `other_seq` where the baseline starts.
    /// This will effectively ignore the first `offset` characters of `other_seq`, as if they were not present.
    /// # Arguments
    /// * `offset` - the number of characters to ignore
    pub fn set_offset(&mut self, offset: usize) {
        self.offset = offset;
    }

    /// This is the main function to extend the WFA with a new symbol.
    /// This function will handle all of the extending and potential edit distance increases as well.
    /// If another edit would exceed the configured cap, the distance freezes at the cap and the DWFA stops voting.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    pub fn update(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Finalized => {
                bail!("Cannot push more bases after finalizing a DWFA");
            },
            DWFALiteState::ExceededEditDistanceLimit => return Ok(self.edit_distance),
            DWFALiteState::Active => {}
        };

        // maximally extend everything along the current diagonals
        let mut maximum_distance = self.extend(baseline_seq, other_seq);
        while maximum_distance < other_seq.len() && !(self.config.allow_early_termination && self.reached_baseline_end(baseline_seq)){
            if self.edit_distance_at_limit() {
                // another edit would pass the cap; freeze and stop contributing
                self.state = DWFALiteState::ExceededEditDistanceLimit;
                return Ok(self.edit_distance);
            }
            // increase the edit distance, re-extension happens automatically
            maximum_distance = self.increase_edit_distance(baseline_seq, other_seq)?;
        }

        // final assertion just to make sure we don't break anything
        assert!(
            maximum_distance == other_seq.len() || 
            (self.config.allow_early_termination && self.maximum_baseline_distance() == baseline_seq.len())
        );
        Ok(self.edit_distance)
    }

    /// This will take the current wavefront and try to extend each one along it's current diagonal.
    /// In most cases, this will not do much work.
    /// This should also be stable such that calling it after finalizing has no impact.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    ///
    /// Returns the maximum reach into `other_seq` after extending every diagonal.
    fn extend(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> usize {
        let mut maximum_distance = 0;

        for (i, d) in self.wavefront.iter_mut().enumerate() {
            // `i` is the index in the wavefront
            // `i // 2` is always the middle diagonal
            // anything less than that has used more baseline symbols
            // anything more than that has used more other symbols
            // anything above or below needs to shift the comparison

            // for this particular wavefront, extend as far as possible
            loop {
                // example at edit distance 1:
                // baseline offset would equal: [d+1, d, d-1] which corresponds to
                // a skipped base in primary, equal usage, and a skipped base in other

                // if we truncate a wavefront (e.g., a bounded size); then `i` below will be too small for every wavefront from the start
                //     that gets truncated; this means we want the below formula to be `*d + self.edit_distance - (i + truncate_prefix)`
                //     where `truncate_prefix` is the number we have cut off the front (in total).
                //     IMPORTANT: this formula will need to get updated everywhere that `baseline_offset` is calculated

                let baseline_offset = *d + self.edit_distance - i;
                let other_offset = *d + self.offset;
                if baseline_offset >= baseline_seq.len() ||
                    other_offset >= other_seq.len() ||
                    (
                        // not equal 
                        baseline_seq[baseline_offset] != other_seq[other_offset] && 
                        // AND baseline is not the configured wildcard
                        self.config.wildcard.is_none_or(|wildcard| baseline_seq[baseline_offset] != wildcard)
                    ) {
                    // if we are past the end of either sequence OR
                    // the sequences are not equal at this position THEN
                    // we are done extending
                    break;
                }

                // we are not done, so add one to this wavefront
                *d += 1;
            }
            maximum_distance = maximum_distance.max(*d);
        }
        self.offset + maximum_distance
    }

    /// This will increase the edit distance for this DWFA and create a new larger wavefront.
    /// This function will automatically call `extend()` to enforce the assumptions.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If the DWFA is already finalized
    fn increase_edit_distance(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Active => {},
            DWFALiteState::Finalized => {
                bail!("Cannot increase edit distance after finalizing a DWFA");
            },
            DWFALiteState::ExceededEditDistanceLimit => {
                bail!("Cannot increase edit distance after exceeding the edit distance limit");
            },
        };

        // first, increase the distance we're at
        self.edit_distance += 1;

        let old_len = self.wavefront.len();
        self.wavefront.resize(old_len + 2, 0);
        // Old slots supply deletion transitions. Reverse traversal reads each
        // old value before mismatch/insertion transitions can overwrite it.
        for i in (0..old_len).rev() {
            let advanced = self.wavefront[i] + 1;
            self.wavefront[i + 1] = self.wavefront[i + 1].max(advanced);
            self.wavefront[i + 2] = self.wavefront[i + 2].max(advanced);
        }

        // re-extend
        Ok(self.extend(baseline_seq, other_seq))
    }

    /// This function signals that base insertion into `other_seq` is completed.
    /// Ordinary mode requires one alignment path to consume both sequences, potentially increasing edit distance.
    /// Early-termination mode requires only the baseline to finish, leaving the other sequence's suffix unpenalized.
    /// Callers must first successfully update with these sequences and a valid offset.
    /// Successful finalization seals the instance against further updates or finalization.
    /// If reaching the required endpoint would exceed the cap, returns successfully in
    /// [`ExceededEditDistanceLimit`](DWFALiteState::ExceededEditDistanceLimit) with the distance frozen at the cap.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If the DWFA is already finalized.
    pub fn finalize(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Finalized => {
                bail!("Cannot finalize a DWFA twice.");
            },
            DWFALiteState::ExceededEditDistanceLimit => return Ok(()),
            DWFALiteState::Active => {},
        };
        while if self.config.allow_early_termination {
            self.maximum_baseline_distance() < baseline_seq.len()
        } else {
            !self.reached_global_end(baseline_seq.len(), other_seq.len())
        } {
            if self.edit_distance_at_limit() {
                self.state = DWFALiteState::ExceededEditDistanceLimit;
                return Ok(());
            }
            self.increase_edit_distance(baseline_seq, other_seq)?;
        }
        self.state = DWFALiteState::Finalized;
        Ok(())
    }

    /// Checks that one path reaches both ends on the terminal diagonal.
    fn reached_global_end(&self, baseline_len: usize, other_len: usize) -> bool {
        let effective_other_len = other_len - self.offset;
        let gap = effective_other_len.abs_diff(baseline_len);
        if gap > self.edit_distance {
            return false;
        }
        let index = if effective_other_len >= baseline_len {
            self.edit_distance + gap
        } else {
            self.edit_distance - gap
        };
        self.wavefront[index] >= effective_other_len
    }

    /// True when the next edit would pass the cap resolved at construction.
    fn edit_distance_at_limit(&self) -> bool {
        self.config.max_edit_distance.is_some_and(|limit| self.edit_distance >= limit)
    }

    /// Helper function that will determine the farthest distance reached into the `baseline_seq` so far.
    pub fn maximum_baseline_distance(&self) -> usize {
        // baseline distance requires some compute
        // the 0-index corresponds to deleting `edit_distance` bases in `baseline`, so it has the largest offset
        // each additional iteration pushes the diagonal closer to inserting bases into `baseline`, so the shift gets progressively smaller
        self.wavefront.iter().enumerate()
            .map(|(i, &d)| d + self.edit_distance - i)
            .max().unwrap()
    }

    /// Helper function that will determine the farthest distance reached into the `other_seq` so far.
    /// After an update, this is the `other_seq` length unless early termination or the edit-distance limit stopped extension.
    pub fn maximum_other_distance(&self) -> usize {
        // other distance is directly tracked in our wavefront
        self.offset + *self.wavefront.iter().max().unwrap()
    }

    /// Returns true if the farther wavefront in the baseline is at the end
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    pub fn reached_baseline_end(&self, baseline_seq: &[u8]) -> bool {
        self.maximum_baseline_distance() == baseline_seq.len()
    }

    /// This will return the set of candidate extensions that do not require increasing the edit distance.
    /// Inactive DWFAs do not nominate extensions.
    /// It also includes how many times that character was counted in the event of multiple possible extension points.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    pub fn get_extension_candidates(&self, baseline_seq: &[u8], other_seq: &[u8]) -> HashMap<u8, usize> {
        let mut ret: HashMap<u8, usize> = Default::default();
        self.for_each_extension_candidate(baseline_seq, other_seq, |symbol| {
            let entry = ret.entry(symbol).or_insert(0);
            *entry += 1;
        });
        ret
    }

    pub(crate) fn fill_extension_candidates(&self, baseline_seq: &[u8], other_seq: &[u8], scratch: &mut CandidateScratch) {
        scratch.reset();
        self.for_each_extension_candidate(baseline_seq, other_seq, |symbol| scratch.record(symbol));
    }

    fn for_each_extension_candidate(&self, baseline_seq: &[u8], other_seq: &[u8], mut record: impl FnMut(u8)) {
        if self.state != DWFALiteState::Active {
            return;
        }
        for (i, &d) in self.wavefront.iter().enumerate() {
            let other_offset = d + self.offset;
            if other_offset == other_seq.len() {
                let offset = d + self.edit_distance - i;
                if offset < baseline_seq.len() {
                    record(baseline_seq[offset]);
                }
            }
        }
    }

    // Getters below
    pub fn edit_distance(&self) -> usize {
        self.edit_distance
    }

    pub fn wavefront(&self) -> &[usize] {
        &self.wavefront
    }

    pub fn max_edit_distance(&self) -> Option<usize> {
        self.config.max_edit_distance
    }

    pub fn state(&self) -> DWFALiteState {
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cdwfa_config::CdwfaConfig;

    fn dwfa_config(wildcard: Option<u8>, allow_early_termination: bool, max_edit_distance: Option<usize>) -> DWFALiteConfig {
        DWFALiteConfigBuilder::default()
            .wildcard(wildcard)
            .allow_early_termination(allow_early_termination)
            .max_edit_distance(max_edit_distance)
            .build()
            .unwrap()
    }

    /// Resolves a consensus cap and stores that absolute maximum on a new DWFA.
    fn dwfa_with_consensus_cap(
        baseline_len: usize,
        hard_max: Option<usize>,
        fraction: Option<f64>,
        floor: Option<usize>,
    ) -> Result<DWFALite, Box<dyn std::error::Error>> {
        let max_edit_distance = CdwfaConfig {
            max_edit_distance: hard_max,
            max_edit_distance_fraction: fraction,
            min_edit_distance: floor,
            ..Default::default()
        }.max_edit_distance_for(baseline_len)?;
        Ok(DWFALite::new(dwfa_config(None, false, max_edit_distance)))
    }

    #[test]
    fn test_finalize_endpoint_modes() {
        let cases: &[(&[u8], &[u8], usize, usize, usize)] = &[
            (b"AC", b"CA", 0, 2, 1),
            (b"AC", b"ACTT", 0, 2, 0),
            (b"ACTT", b"AC", 0, 2, 2),
            (b"AC", b"GGCA", 2, 2, 1),
            (b"AC", b"GGACTT", 2, 2, 0),
        ];
        for &(baseline, other, offset, global_score, early_score) in cases {
            for early in [false, true] {
                for incremental in [false, true] {
                    let mut dwfa = DWFALite::new(dwfa_config(None, early, None));
                    dwfa.set_offset(offset);
                    if incremental {
                        for length in offset..=other.len() {
                            dwfa.update(baseline, &other[..length]).unwrap();
                        }
                    } else {
                        dwfa.update(baseline, other).unwrap();
                    }
                    dwfa.finalize(baseline, other).unwrap();
                    assert_eq!(dwfa.state(), DWFALiteState::Finalized);
                    assert_eq!(
                        dwfa.edit_distance(), if early { early_score } else { global_score },
                        "baseline={baseline:?}, other={other:?}, offset={offset}, early={early}, incremental={incremental}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_finalize_matches_pairwise_oracle() {
        let sequences: Vec<Vec<u8>> = (0..=4)
            .flat_map(|length| {
                (0..1usize << length).map(move |bits| {
                    (0..length).map(|i| b"AC"[(bits >> i) & 1]).collect()
                })
            })
            .collect();
        for baseline in &sequences {
            for suffix in &sequences {
                for offset in [0, 2] {
                    let mut other = vec![b'G'; offset];
                    other.extend_from_slice(suffix);
                    for early in [false, true] {
                        let expected = crate::sequence_alignment::wfa_ed_config(
                            &other[offset..], baseline, !early, None
                        );
                        for incremental in [false, true] {
                            let mut dwfa = DWFALite::new(dwfa_config(None, early, None));
                            dwfa.set_offset(offset);
                            if incremental {
                                for length in offset..=other.len() {
                                    dwfa.update(baseline, &other[..length]).unwrap();
                                }
                            } else {
                                dwfa.update(baseline, &other).unwrap();
                            }
                            dwfa.finalize(baseline, &other).unwrap();
                            assert_eq!(
                                dwfa.edit_distance(), expected,
                                "baseline={baseline:?}, other={other:?}, offset={offset}, early={early}, incremental={incremental}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_finalize_clone_isolation() {
        let mut live = DWFALite::default();
        live.update(b"AC", b"A").unwrap();
        let mut finalized = live.clone();
        finalized.finalize(b"AC", b"A").unwrap();
        assert_eq!(finalized.edit_distance(), 1);

        live.update(b"AC", b"AC").unwrap();
        live.finalize(b"AC", b"AC").unwrap();
        assert_eq!(live.edit_distance(), 0);

        let mut sealed_clone = finalized.clone();
        assert!(sealed_clone.finalize(b"AC", b"A").is_err());
        assert_eq!(sealed_clone, finalized);
        assert!(sealed_clone.update(b"AC", b"AC").is_err());
        assert_eq!(sealed_clone, finalized);
    }

    #[test]
    fn test_finalize_global_endpoint() {
        let mut dwfa = DWFALite::default();
        assert_eq!(dwfa.update(b"AC", b"CA").unwrap(), 1);
        dwfa.finalize(b"AC", b"CA").unwrap();
        assert_eq!(dwfa.edit_distance(), 2);
    }

    #[test]
    fn test_finalize_seals_alignment() {
        for early in [false, true] {
            for (baseline, expected) in [(b"AC".as_slice(), 0), (b"ACTT".as_slice(), 2)] {
                let mut dwfa = DWFALite::new(dwfa_config(None, early, None));
                assert_eq!(dwfa.update(baseline, b"AC").unwrap(), 0);
                dwfa.finalize(baseline, b"AC").unwrap();
                assert_eq!(dwfa.edit_distance(), expected);
                assert_eq!(dwfa.state(), DWFALiteState::Finalized);
                let snapshot = dwfa.clone();
                assert!(dwfa.finalize(baseline, b"AC").is_err());
                assert_eq!(dwfa, snapshot);
                assert!(dwfa.update(baseline, b"ACA").is_err());
                assert_eq!(dwfa, snapshot);
            }
        }
    }

    #[test]
    fn test_incremental_reach_with_offset_and_early_termination() {
        let baseline = b"AC";
        let other = b"GGACTT";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.set_offset(2);
        for length in 2..=other.len() {
            assert_eq!(dwfa.update(baseline, &other[..length]).unwrap(), 0);
            assert_eq!(dwfa.maximum_other_distance(), length.min(4));
            assert_eq!(dwfa.maximum_baseline_distance(), (length - 2).min(2));
        }
        dwfa.finalize(baseline, other).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_new() {
        let dwfa = DWFALite::default();
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa.wavefront(), &[0]);
    }

    #[test]
    fn test_exact_match() {
        let sequence = b"ACGTACGTACGT";
        let mut other_seq = vec![];
        let mut dwfa = DWFALite::default();
        for &c in sequence.iter() {
            other_seq.push(c);
            assert_eq!(dwfa.update(sequence, &other_seq).unwrap(), 0);
        }
    }

    #[test]
    fn test_simple_mismatch() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACCTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_candidate_counts_across_edit_layers() {
        let baseline = [0, 255, 0];
        let other = [255, 0, 0];
        let expected: [(usize, &[(u8, usize)]); 3] = [
            (1, &[(0, 2), (255, 1)]),
            (1, &[(255, 1)]),
            (2, &[(0, 1), (255, 1)])
        ];
        let mut dwfa = DWFALite::default();
        let mut scratch = CandidateScratch::default();
        for (index, (score, counts)) in expected.iter().enumerate() {
            let prefix = &other[..=index];
            assert_eq!(dwfa.update(&baseline, prefix).unwrap(), *score);
            let expected_counts: HashMap<u8, usize> = counts.iter().copied().collect();
            assert_eq!(dwfa.get_extension_candidates(&baseline, prefix), expected_counts);
            dwfa.fill_extension_candidates(&baseline, prefix, &mut scratch);
            assert_eq!(scratch.ordered_counts().collect::<HashMap<_, _>>(), expected_counts);
        }
        dwfa.finalize(&baseline, &other).unwrap();
        assert_eq!(dwfa.edit_distance(), 2);
    }

    #[test]
    fn test_simple_insertion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACIGTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_simple_deletion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_complex_001() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACTACGCACGGGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 4);
    }

    #[test]
    fn test_complex_002() {
        //modified_seq has 2 separate deletions, 1 2bp insertion, and 1 mismatch
        let sequence     = b"AACGGATCAAGCTTACCAGTATTTACGT";
        let alt_sequence = b"AACGGACAAAAGCTTACCTGTATTACGT";
        
        let mut dwfa = DWFALite::default();
        /*
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        */
        // test a shortcut
        dwfa.update(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), 5);
    }

    #[test]
    fn test_big_insertion() {
        // one big insertion in the middle
        let sequence     = b"AACGGATTTTACGT";
        let alt_sequence = b"AACGGATAAAAGCTTACCTGTTTTACGT";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), alt_sequence.len() - sequence.len());
    }

    #[test]
    fn test_big_deletion() {
        // one big deletion in the middle
        let sequence     = b"ATTTTTTTTTTAAAAAAAAAA";
        let alt_sequence = b"AAAAAAAAAAA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_required_finalize() {
        // the ALT is a lot smaller than the baseline, so we need to run a finalize to get the full ED
        let sequence     = b"ATTTTTTTTTTA";
        let alt_sequence = b"AA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }

        // here it has only compared AT to AA, so ED=1
        assert_eq!(dwfa.edit_distance(), 1);

        // after finalizing, it will do end-to-end comparison
        dwfa.finalize(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_cloning() {
        let sequence     = b"AAAAAAA";
        let alt_sequence = b"AAACAAA";

        let mut dwfa = DWFALite::default();
        let mut dwfa2 = dwfa.clone();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &sequence[..(l+1)]).unwrap();
            dwfa2.update(sequence, &alt_sequence[..(l+1)]).unwrap();

            if sequence[l] == alt_sequence[l] {
                // same sequence still
                assert_eq!(dwfa, dwfa2);
            } else {
                // should have different sequences now
                assert_ne!(dwfa, dwfa2);

                // re-clone the first one
                dwfa2 = dwfa.clone();
            }
        }

        // in the end, both should exactly match due to cloning
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa2.edit_distance(), 0);
    }

    #[test]
    fn test_wildcards_001() {
        //modified_seq has several wildcards, but otherwise an exact match
        let consensus= b"AACGGATCAAGCTTACCAGTATTTACGT";
        let baseline = b"*ACGGATCAA**TTACCA*TATTTACG*";
        
        let mut dwfa = DWFALite::new(dwfa_config(Some(b'*'), false, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_wildcards_002() {
        //modified_seq has several wildcards, and we added 1 del, 1 SNV, and 1 ins
        let consensus= b"AACGGATCAAGCTTACCAGTATTTACGT";
        let baseline = b"*ACGATCAA**TATACCA*TATCTACG*";
        
        let mut dwfa = DWFALite::new(dwfa_config(Some(b'*'), false, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 3);
    }

    #[test]
    fn test_early_termination_001() {
        let consensus = b"ACGTACGT";
        let baseline = b"ACGT";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_big_early_termination() {
        let c1 = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACACGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCATAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGTTCGGGGTGGAAGCGGAAGGGCTTCTCCCAGACGGCCTCATCCTTCAGCACCGATGACAGGTTGGTGATGAGTGTCGTTCCCTGGGCAGGAGATGCAGGGTGAGAGTGGGGACTGGACTCTAGGATGCTGGGACCCCTGCCACCAAACACACGGGGGACACACACTGCCTGGCACACAGCTGGACTCTGTCAACTAGTCCTGCGCCCGAGAAGCTCCACAGTACCCTCTCCGACCCCACAGCAGGGCGCAGTCACACCTCTCAGAGGCACCCACACTGCCCCCTCTCCCTGCAGGCGCTGGGTCCTCCAACATTCTGGCAGGTCCTGGTTTGTCTCCCCACTAGACGGGGGCTCTGGATGGACAGGCCAGCCCTGCCTATACTCTGGACCCCCCACCCAAGTGGGGACAGTCAGTGTGGTGGCATTGAGGACTAGGTGGCCAGGGTTCCTAGAGTGGGCCCACCTGGCAGTAGCCATGCTGGGGCTATCACCAGGGGCTGGTGCTGAGCTGGGGTGAGGAGGGCGCCAGGCCTACCTTAGGGATGCGGAAGCCCTGTACTTCGATGTCACGGGATGTCATATGGGTCACACCCAGGGGGACGATGTCCCCAAAGCGCTGCACCTCATGAATCACGGCAGTGGTGTAGGGCATGTGAGCCTGGTCACCCATCTCTGGTCGCCGCACCTGCCCTATCACGTCGTCGATCTCCTGTTGGACACGGACTGGACAGACATGCGTCCCCACAATGGGTCAGCACCCAGGGGACACTCTCCTTCCTCCTGTGTTGGAGGAAGTTAGGCTTACAGGAGCCTGGCCACGCCTGTGCTGGAAGCCCCGGGTGTCCCAGCTAAGCCCAGGGGCCCCCAGCTGTACCCTTCCTCCCTCAGTCCCTGCCTTGGGCCCCAGCTGGGCTCACGCTGCACATCCAGGTGTAGGATCATGAGCAGGAGGCCCCAGGCCAGCGTGGTCAAGGTGGTCACCATCCCGGCAAGGAACAGGTTACCCACCACTATGCGCAGGTTCTCATCATTGAAGCTGCTCTCAGGGCTCCCCTTGGCCTGAGCAGGGCCGAGAGGATACTCAGGGGATAGAACGGGGTAGCCCCCAAATGACCTCCAATTCTGCACCTGTCAGCCCAGATGCGGCTCGCCGGGTGATGCACTGGTCCAACCTTTTGCCCAGCCTCCCCTCATTCCTCCTGGGACGTTCAACCCACCACCCTTGCCCCCCACCGTGGCAGCCACTCTCACCTTCTCCTTCTTTGCCAGGAAGGCCTCAGTCAGGTCTCGGGGTGGCTGGGCTGGGTCCCAGGTCATCCTGTGCTCAGTTAGCAGCTCATCCAGCTGGGTCAGGAAAGCCTTTTGGAAGCGTAGGACCTTGCCAGCCAGCGCTGGGATGTGCGGGAGGACGGGGACAGCATTCAGCACCTACACCAGACAGAACCGGGTCTCAATCCTTCCTGTGCTCTGCGTTCATCTGGACCAGTCTCAGGCCCCAGCCATCTCCAGGAAGACCCAGGGCCTGCCTGTCCTTACCACTGACCTCACCAAGTCCCTCCCCAAGTGCCAGCCTCCACCCTCTCTCTCCTTGCCCAGAGGAGAAACCTAAAATCGAAATCTCCAACGTGGACGGGGGTACAGAGTCCTTGGCCTCTCCTGGTGCCCCCTGACCCGGGCACACCTCTCCCACGACCATGTCTGAGATGTCCCCTCCTCCTCCAGGCCCTTCTTACAGTGGGGTCTCCTGGAATGTCCTTTCCCAAACCCATCTACGCAAATCCTGCCCTTCGGAGGCCCCAGTCCAGCCCCGGCACCTCTCAGGAGCTCGCCCTGCAAAGACCCTTGCTCCGCACCTCGCGCAGGAAGCCCGACTCCTCCTTCGATCCCTCCCTGAGCTAGGTCCAGCAGCCTGAGGAAGCGAGGGTCGTCGTACTCGAAGCGGCGCCCGCAGGTGAGGGAGGCGATCACGTTGCTCACGGCTTTGTCCAAGAGACCGTTGGGGCGAAAGGGGCGTCCTGGGGGTGGGAGATGCGGGTAAGGGGTCGCCTTCTCCGTCCCCCGCCTTCCCAGTTCCCGCTGTGTGCCCTTCTGCCCATCACCCACCGGCTTGGTCGGCGAAGGCGGCACAAAGGCAGGCGGCCTCCTCGGTCACCCACTGCTCCAGCGACTTCTTGCCCAGGCCCAAGTTGCGCAAGGTGGACACGGAGAAGCGCCTCTGCTCGCGCCACGCGGGCCCATAGCGCGACAGGATCACCCCTGTGGGCGGGACGGACACGTGGGCGTTGCCATGAAGGCCTTGGCCCCACCCTCCGCCACCCACTCCAACCCTGGCGCTCCACAAGGTCTCCCGCAGTCCCTAGCCCGGTCCAGCTGGGCACAGGGCCCACTCTTTGCTCACCCACATTGCTCCCCTGCCTGGGGCGGGGTTTGGCCCCACCTCGTCTCTGCCCACCCTGACCACCTTTCCACTCAAGGAAGATCCCGCCCGTCCCGCCCACACTGAGCCCGCAGCATAGGCGCGGTCCCCGCCACCGCCACTTCGACGCATCAGCCTCGCCCACCGGGCTTCTGGCGGGTCTGGGCAGTAGCCCCGCCCCCTCCCAGCCCACAGACTCGCACCTCCCCCGTGCAGGTGGTTTCCTGGCCCACTGTCCTCAGCCCACTCGCTGGCCTTTATCTCTGTTTCACGTCCAGGACCCCACGCCCTGTCGGCGCTGCTTGGGCTACGGTCACTGTCCACCCGGGGCCCACGGAAACGCGGTCTCTGTCCCCCACCGCCGCTTGCCTTGGGAACGCGGCCCGAAGCCCAGGACCTGGTAGATGGGCGCAGGCGGGCGGTCGGCCGTGTCCTCGCCGCGGGTCACCATCGCCTCGCGCACGGCCGCCAGCCCATTGAGCACGACCACCGGCGTCCAGGCCAGCTGCAGGCTGAACACGTCCCCGAAGCGGCGCCGCAACTGCAGAGGGAGGGTCAGGGCCTCTTGTCAAGCCAGGATCCCCCCAGACTACAGGTCCTAGTCCTATTTGAACCTTGGACGACCCCCGGGGCTACCAGGAGTGAGCAGGTGGAAGGAGGAGACCCAGCCTCCTGATCCTGGGGCGGGGGTGGGGGTCACACCTTCTGTGATGGAGGAACTCAGTTTGGATGCGTCACCCAGGTATGACCTTGCAAGAGTCACCAAAATTGCCGAGAGGCCCCAGTTAGCATCCCATTCCCAGATGATGGTCCATGCCGGTGAGCAGTGAGGCCCGAGGACCCACAGTGCAAAAGGTTTGAACCGGGTCACTGCACCCCCTTCATCCTCGATTTCGTGATTTAAACGGCACTCAGGACTAACTCATCTTCCATTCCCAAGGCCTTTCCTTCTGGTGTCAGCAGAAGGGACTTTGTACTCCATAACATATGTTGCCCAATGGGCTTGCATGCCCACTGCCAAGTCCAGCTCCACCTCCAGGCCCTTGCCCTACTCTTCCTTGGCCTTTGGAAAATCCAGTCCTTCATGCCATGTATAAATGTCCTTCCCCAGGACGTCCCCCAAACCTGCTTCCCCTTCTCAGCCTGGCTTCTGATCCAGCCTGTGGTTTAACCCACCACCCATGTTTGCTGGTGGTGGGGCATCCTCAGGACCTCTGCCGCCCTCCAGGACCTCCTCCCTCACCTGGTCGAAGCAGTATGGTGTGTTCTGGAAGTCCACATGCAGCAAGGTTGCCCAGCCCGGGCAGTGGCAGGGGACCTGGCGGGTAGCGTGCAGCCCAGCGTTGGTGCCGGTGCATCAGGTCCACCAGGAGCAGGAAGATGGCCACTATCATGGCCAGGGGCACCAGTGCTTCTAGCCCCATGGCTGCCTCACTACCAACTGGGCTCCTCTGGACACACCTGGCACCCCCACCCCACCAGGCACAGAGGACCAGGCAGGACACTCTCAGCACACCGAGCGCGTGACCCTTCCCTTATAAAGGGAGCTGATGATGGCCTTCGCCCTCTGCTGTGAGTGAACCTGCTGTGTTGACTGTGCTGCCAGTGGCAGAGTCAGGCCAGGGTGGGTATGGGCTGCTCCAGAGGTCCTTGCCGCTGCTTCCTGCTCCAGGCCCTTACCCAGGGTAGGGTGGTAGAAAGGCCTGGTCGGAGAAGTCACCCCCTCTCCCCACTCCAAGCTCCCCAAGCCCACACAGGCTTCTGGGATAACCAGGGTCTCAGTGGACCCGGCCATCCACCTCCCAGCTAGGCTCATACACCGTAATGTAGTCACAACCCCTCCTCCAGAACATGGCCTTGCCCTTTCCCTACCCCCACCTGCCCACTCCAGAGTGACCTTCAGCACCCTTATCTGTCACTGGCACTTACCTGGGGCCTTAGAGCTCCTGATGATGAGTGGCATCATGGGCCTGGTCCCTTCACTTCACCTTGCACTCTTGACATGCACAGACGCTATGCACACACCTGATGGTGCACAGATCTCTTGTCCACTCCCAGACACTTGTCCACTTGTTCACACTTGCAGGGACACGATTACACATGCAGAAAATCACCCACACAAAGACAATATTCACACATACACAGACTCACACTGACACTCAGGGCACACATTCTCTCTCACACACACCAGTCACACACACATACAGACCCGGCACCAAGTACCCCACTTCCCAGCCATGCCCAAGGTTTCCTGGATGGGACCTCTCCTGTCCAGAGGCTGCTCCCAGTGAGCCTCAAAGCTGTCACGTGGATCCCAGCTCAGCCCACATTCTGGGCTCTGGCCGGGCCATGGCTTCCTGTTTGCAACAGGGCTGTTCCCAGAGCTCCCAGTTGGTAGCCTGAAGGCCCTTGCCCCAGCCTGTGACAGCATCCTCCAGGGCTGCCTGAGGGTCGTCATTCTCCACTGCTTCCTGGCCTCCATGTTTCTGATTAGAAATCTGGTGGAAACATTATGGAGGATCCTTTATTTAGGATATGTTGCTTTTTTATTTTTATTTTTTCTTTAGACAGGGTCTCACTCTGTTGCCCGGGCCGGAGTGCAGTGGCAGGATCACGGCTCACTGCAATCTCAACATCAAGTGGACCTCCTGCCTCCCAAGTAGCTGGGACTACAGGCACCACCGAGCCCAAATAATTTTTTTTTTGAGACGGAGTTTTGCTCTGTCGCCCAGGTGGGAGTGCAATGATGCGATCTCGGCTCACTGCAACCTCCACCTCCAGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAGGTGCCCACCACCATGCCTGGCTGATTTTTTGTA";
        let seq_23 = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCGTAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGT";

        // iterate, making sure everything is fine even as we go well beyond seq_23
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        for i in 0..c1.len() {
            dwfa.update(seq_23.as_bytes(), c1[0..(i+1)].as_bytes()).unwrap();
            assert!(dwfa.edit_distance() <= 2);
        }
        assert_eq!(dwfa.edit_distance(), 2);

        // make sure finalize does not break it
        dwfa.finalize(seq_23.as_bytes(), c1.as_bytes()).unwrap();
        assert_eq!(dwfa.edit_distance(), 2);
        
    }

    #[test]
    fn test_offsets() {
        // this sequence
        let consensus = b"ACGTACGT";
        let baseline =    b"GTACGT";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.set_offset(2);
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_exact_match_under_edit_distance_cap() {
        let sequence = b"ACGTACGT";
        let mut dwfa = dwfa_with_consensus_cap(sequence.len(), Some(2), None, None).unwrap();
        let other = &sequence[..4];
        assert_eq!(dwfa.update(sequence, other).unwrap(), 0);
        assert_eq!(dwfa.get_extension_candidates(sequence, other).get(&b'A'), Some(&1));
        assert_eq!(dwfa.state(), DWFALiteState::Active);
    }

    #[test]
    fn test_finalize_global_endpoint_under_edit_distance_cap() {
        let baseline = b"AC";
        for other in [b"CA".as_slice(), b"GGCA".as_slice()] {
            let offset = other.len() - 2;
            for early in [false, true] {
                for cap in [1, 2] {
                    for incremental in [false, true] {
                        let mut dwfa = DWFALite::new(dwfa_config(None, early, Some(cap)));
                        dwfa.set_offset(offset);
                        if incremental {
                            for length in offset..=other.len() {
                                dwfa.update(baseline, &other[..length]).unwrap();
                            }
                        } else {
                            dwfa.update(baseline, other).unwrap();
                        }
                        assert_eq!(dwfa.edit_distance(), 1);
                        assert_eq!(dwfa.state(), DWFALiteState::Active);
                        // Different paths reach each end at score 1, but no single path reaches both.
                        assert_eq!(dwfa.maximum_baseline_distance(), baseline.len());
                        assert_eq!(dwfa.maximum_other_distance(), other.len());
                        assert!(!dwfa.get_extension_candidates(baseline, other).is_empty());

                        dwfa.finalize(baseline, other).unwrap();
                        let excluded = !early && cap == 1;
                        assert_eq!(dwfa.edit_distance(), if early { 1 } else { cap });
                        assert_eq!(
                            dwfa.state(),
                            if excluded { DWFALiteState::ExceededEditDistanceLimit } else { DWFALiteState::Finalized }
                        );
                        assert!(dwfa.get_extension_candidates(baseline, other).is_empty());
                        let mut scratch = CandidateScratch::default();
                        scratch.record(b'X');
                        dwfa.fill_extension_candidates(baseline, other, &mut scratch);
                        assert!(scratch.ordered_counts().next().is_none());

                        if excluded {
                            let snapshot = dwfa.clone();
                            let mut longer = other.to_vec();
                            longer.push(b'C');
                            assert_eq!(dwfa.update(baseline, &longer).unwrap(), cap);
                            dwfa.finalize(baseline, &longer).unwrap();
                            assert_eq!(dwfa, snapshot);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_finalize_cap_while_consuming_baseline_suffix() {
        for early in [false, true] {
            for cap in [1, 2] {
                let mut dwfa = DWFALite::new(dwfa_config(None, early, Some(cap)));
                assert_eq!(dwfa.update(b"ACTT", b"AC").unwrap(), 0);
                dwfa.finalize(b"ACTT", b"AC").unwrap();
                assert_eq!(dwfa.edit_distance(), cap);
                assert_eq!(
                    dwfa.state(),
                    if cap == 1 { DWFALiteState::ExceededEditDistanceLimit } else { DWFALiteState::Finalized }
                );
                assert!(dwfa.get_extension_candidates(b"ACTT", b"AC").is_empty());
                let mut scratch = CandidateScratch::default();
                scratch.record(b'T');
                dwfa.fill_extension_candidates(b"ACTT", b"AC", &mut scratch);
                assert!(scratch.ordered_counts().next().is_none());
            }
        }
    }

    #[test]
    fn test_zero_edit_distance_cap_endpoint_and_early_termination() {
        for early in [false, true] {
            let mut dwfa = DWFALite::new(dwfa_config(None, early, Some(0)));
            assert_eq!(dwfa.update(b"AC", b"A").unwrap(), 0);
            assert_eq!(dwfa.state(), DWFALiteState::Active);
            assert_eq!(dwfa.get_extension_candidates(b"AC", b"A").get(&b'C'), Some(&1));
            assert_eq!(dwfa.update(b"AC", b"AC").unwrap(), 0);
            assert_eq!(dwfa.state(), DWFALiteState::Active);
            let mut finalized = dwfa.clone();
            finalized.finalize(b"AC", b"AC").unwrap();
            assert_eq!(finalized.state(), DWFALiteState::Finalized);
            assert_eq!(finalized.edit_distance(), 0);

            assert_eq!(dwfa.update(b"AC", b"ACTT").unwrap(), 0);
            assert_eq!(
                dwfa.state(),
                if early { DWFALiteState::Active } else { DWFALiteState::ExceededEditDistanceLimit }
            );
            dwfa.finalize(b"AC", b"ACTT").unwrap();
            assert_eq!(
                dwfa.state(),
                if early { DWFALiteState::Finalized } else { DWFALiteState::ExceededEditDistanceLimit }
            );
            assert_eq!(dwfa.edit_distance(), 0);
        }
    }

    #[test]
    fn test_hard_edit_distance_cap() {
        let baseline = b"ACGTACGTACGT";
        let other = b"TTTTTTTTTTTT";
        let mut dwfa = dwfa_with_consensus_cap(baseline.len(), Some(3), None, None).unwrap();
        assert_eq!(dwfa.update(baseline, other).unwrap(), 3);
        assert_eq!(dwfa.state(), DWFALiteState::ExceededEditDistanceLimit);
        assert!(dwfa.get_extension_candidates(baseline, other).is_empty());
        let mut scratch = CandidateScratch::default();
        scratch.record(b'A');
        dwfa.fill_extension_candidates(baseline, other, &mut scratch);
        assert!(scratch.ordered_counts().next().is_none());
        let snapshot = dwfa.clone();

        let mut longer = other.to_vec();
        longer.push(b'A');
        assert_eq!(dwfa.update(baseline, &longer).unwrap(), 3);
        dwfa.finalize(baseline, &longer).unwrap();
        assert_eq!(dwfa.edit_distance(), 3);
        assert!(dwfa.get_extension_candidates(baseline, &longer).is_empty());
        dwfa.fill_extension_candidates(baseline, &longer, &mut scratch);
        assert!(scratch.ordered_counts().next().is_none());
        assert_eq!(dwfa, snapshot);
    }

    #[test]
    fn test_fractional_edit_distance_cap() {
        let baseline = vec![b'A'; 100];
        let mut other = vec![b'T'; 6];
        other.extend(std::iter::repeat(b'A').take(94));
        let mut dwfa = dwfa_with_consensus_cap(baseline.len(), None, Some(0.05), None).unwrap();
        assert_eq!(dwfa.max_edit_distance(), Some(5));
        assert_eq!(dwfa.update(&baseline, &other).unwrap(), 5);
        assert_eq!(dwfa.state(), DWFALiteState::ExceededEditDistanceLimit);
        assert!(dwfa.get_extension_candidates(&baseline, &other).is_empty());
    }

    #[test]
    fn test_edit_distance_cap_uses_tighter_limit() {
        let tighter_hard = dwfa_with_consensus_cap(100, Some(3), Some(0.05), None).unwrap();
        assert_eq!(tighter_hard.max_edit_distance(), Some(3));
        let tighter_fraction = dwfa_with_consensus_cap(100, Some(10), Some(0.05), None).unwrap();
        assert_eq!(tighter_fraction.max_edit_distance(), Some(5));
        let unlimited = dwfa_with_consensus_cap(100, None, None, None).unwrap();
        assert_eq!(unlimited.max_edit_distance(), None);
        assert!(dwfa_with_consensus_cap(100, None, Some(-0.1), None).is_err());
        assert!(dwfa_with_consensus_cap(100, None, Some(f64::NAN), None).is_err());
    }

    #[test]
    fn test_edit_distance_floor_raises_derived_cap() {
        // floor(0.05 * 19) == 0, so the floor is the resolved cap
        let raised = dwfa_with_consensus_cap(19, None, Some(0.05), Some(2)).unwrap();
        assert_eq!(raised.max_edit_distance(), Some(2));
        // a derived cap already above the floor stays at the derived value
        let unchanged = dwfa_with_consensus_cap(100, None, Some(0.05), Some(2)).unwrap();
        assert_eq!(unchanged.max_edit_distance(), Some(5));
        // a floor with no maximum does not create a cap
        let unlimited = dwfa_with_consensus_cap(19, None, None, Some(2)).unwrap();
        assert_eq!(unlimited.max_edit_distance(), None);
    }
}
