
/*!
This module provides access to the ConsensusDWFA, which generates the single best consensus for a set of sequences.

# Example usage
```rust
use waffle_con::consensus::ConsensusDWFA;

let sequences = [
    b"ACGT".to_vec(),
    b"ACCGT".to_vec(), // this should be the consensus
    b"ACCCGT".to_vec()
];

// add all the sequences
let mut cdwfa: ConsensusDWFA = Default::default();
for s in sequences.iter() {
    cdwfa.add_sequence(s).unwrap();
}

// run consensus and check the results
let consensuses = cdwfa.consensus().unwrap();
assert_eq!(consensuses.len(), 1);
assert_eq!(consensuses[0].sequence(), sequences[1]);
assert_eq!(consensuses[0].scores(), &[1, 0, 1]);
assert_eq!(consensuses[0].assignments(), None);
```
*/

use log::{debug, trace};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::{bail, SimpleError};
use std::cmp::{Ordering, Reverse};

use crate::cdwfa_config::{CdwfaConfig, ConsensusCost};
use crate::dynamic_wfa::{DWFALite, DWFALiteState, ExtensionVotes, VoteScratch};
use crate::global_scores::{Columns, Patterns};
use crate::pqueue_tracker::PQueueTracker;
use crate::priority_queue::{CachePriority, PriorityQueue};

/// Lower cost, then longer consensus, then an earlier tie-breaker.
type NodePriority = (Reverse<usize>, usize, Reverse<u64>);

/// Whether a read in a single consensus was used or halted at the edit-distance cap.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConsensusAssignment {
    /// The read was tracked through the consensus.
    Included,
    /// The read stopped contributing because it reached the edit-distance cap.
    EditDistanceLimit,
}

/// Contains a final consensus result
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Consensus {
    /// The generated consensus
    sequence: Vec<u8>,
    /// The consensus scoring model
    consensus_cost: ConsensusCost,
    /// Vector of the scores from the consensus to each sequence
    scores: Vec<usize>,
    /// Per-read assignment. `None` when every score belongs to an included read.
    /// `Some` when any read was halted at the edit-distance cap; the vector matches `scores`.
    assignments: Option<Vec<ConsensusAssignment>>,
}

impl Consensus {
    /// Constructor
    /// # Errors
    /// * if `assignments` is `Some` and its length does not match `scores`
    pub fn new(
        sequence: Vec<u8>,
        consensus_cost: ConsensusCost,
        scores: Vec<usize>,
        assignments: Option<Vec<ConsensusAssignment>>,
    ) -> Result<Consensus, SimpleError> {
        if let Some(assignments) = &assignments {
            if assignments.len() != scores.len() {
                bail!("assignments and scores must have the same length");
            }
        }
        Ok(Consensus {
            sequence,
            consensus_cost,
            scores,
            assignments,
        })
    }

    // Getters
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn consensus_cost(&self) -> ConsensusCost {
        self.consensus_cost
    }

    pub fn scores(&self) -> &[usize] {
        &self.scores
    }

    /// Per-read assignment. `None` when every score belongs to an included read.
    pub fn assignments(&self) -> Option<&[ConsensusAssignment]> {
        self.assignments.as_deref()
    }
}

/// Core utility that will generate a consensus sequence using the DWFA approach.
/// For now, it will assume that all sequences are full length representations
#[derive(Debug, Default)]
pub struct ConsensusDWFA<'a>  {
    /// Contains all the sequences that have been added to this consensus so far.
    sequences: Vec<&'a [u8]>,
    /// Last offset into the consensus that the sequence is allowed to start, this may get adjusted at run-time. If None, then assume start.
    offsets: Vec<Option<usize>>,
    /// The config for this consensus run
    config: CdwfaConfig,
    /// The alphabet we will use for consensus building
    alphabet: HashSet<u8>
}

impl<'a> ConsensusDWFA<'a> {
    /// Creates a new instance of ConsensusDWFA and performs sanity checks.
    /// # Arguments
    /// * `config` - the type of consensus score we want to use
    /// # Errors
    /// * None so far
    pub fn with_config(config: CdwfaConfig) -> Result<ConsensusDWFA<'a>, Box<dyn std::error::Error>> {
        Ok(ConsensusDWFA {
            sequences: vec![],
            offsets: vec![],
            config,
            alphabet: Default::default()
        })
    }

    /// Adds a new sequence to the list without an offset.
    /// # Arguments
    /// * `sequence` - the new sequence to add
    /// # Errors
    /// * None so far
    pub fn add_sequence(&mut self, sequence: &'a [u8]) -> Result<(), Box<dyn std::error::Error>> {
        self.add_sequence_offset(sequence, None)
    }

    /// Adds a new sequence to the list 
    /// # Arguments
    /// * `sequence` - the new sequence to add
    /// * `last_offset` - an optional offset into the consensus that is the last allowed offset for the sequence to start
    /// # Errors
    /// * None so far
    pub fn add_sequence_offset(&mut self, sequence: &'a [u8], last_offset: Option<usize>) -> Result<(), Box<dyn std::error::Error>> {
        // add ever character to our alphabet
        self.alphabet.extend(sequence.iter().cloned());
        
        // make sure we didn't add the wildcard by accident
        // note: it's easier to just remove it once here than do something fancy with the above iteration
        if let Some(wc) = self.config.wildcard {
            self.alphabet.remove(&wc);
        }

        // now add the sequence
        self.sequences.push(sequence);
        self.offsets.push(last_offset);
        Ok(())
    }

    /// The core function that gets called after adding all the sequences we care about
    /// # Errors
    /// * if the DWFA throws errors
    pub fn consensus(&self) -> Result<Vec<Consensus>, Box<dyn std::error::Error>> {
        self.consensus_search(
            (),
            |_, _| {},
            |node, results| {
                // Finalize a clone so the original wavefronts remain available for extension.
                let mut finalized_node = (*node).clone();
                finalized_node.finalize(&self.sequences)?;
                let cost = self.config.consensus_cost;
                let score = finalized_node.total_cost(cost);
                results.retain(score, self.config.max_return_size, || {
                    let scores = finalized_node.costs(cost);
                    let assignments = finalized_node.assignments();
                    Consensus::new(finalized_node.consensus, cost, scores, assignments)
                })
            },
        )
    }

    /// Like [`Self::consensus`], but scores terminals with incremental global edit distances.
    ///
    /// This updates a score column for every read on every extension, even if that branch
    /// is later pruned. Prefer it when finishing wavefronts at many possible endpoints
    /// dominates the search; otherwise [`Self::consensus`] avoids that per-extension work.
    ///
    /// The DWFA still decides which prefixes to extend and score. The incremental scorer
    /// supports L1 scoring without offsets, a wildcard, early termination, or edit-distance
    /// cap settings. Other configurations use [`Self::consensus`] directly.
    ///
    /// # Errors
    /// * if the DWFA throws errors
    pub fn consensus_with_global_scores(&self) -> Result<Vec<Consensus>, Box<dyn std::error::Error>> {
        let config = &self.config;
        if config.consensus_cost != ConsensusCost::L1Distance
            || config.wildcard.is_some()
            || config.allow_early_termination
            || config.max_edit_distance.is_some()
            || config.max_edit_distance_fraction.is_some()
            || config.min_edit_distance.is_some()
            || self.offsets.iter().any(Option::is_some)
        {
            return self.consensus();
        }

        let patterns = Patterns::new(&self.sequences, &self.alphabet);
        self.consensus_search(
            Columns::root(&patterns),
            |scores, symbol| scores.push(&patterns, symbol),
            |node, results| {
                let scores = node.terminal_scores.distances();
                results.retain(scores.iter().sum(), config.max_return_size, || {
                    Consensus::new(
                        node.consensus().to_vec(),
                        ConsensusCost::L1Distance,
                        scores.to_vec(),
                        None,
                    )
                })
            },
        )
    }

    // Specialize the same search for both scorers. The legacy state and append callback
    // are zero-sized; no per-node score storage or dynamic dispatch is needed on that path.
    fn consensus_search<S: Clone>(
        &self,
        terminal_scores: S,
        append_score_column: impl Fn(&mut S, u8),
        consider_terminal: impl Fn(
            &ConsensusNode<S>,
            &mut TerminalResults,
        ) -> Result<(), Box<dyn std::error::Error>>,
    ) -> Result<Vec<Consensus>, Box<dyn std::error::Error>> {
        // initialize everything
        let mut results = TerminalResults { maximum_error: usize::MAX, retained: vec![] };
        let mut nodes_explored: usize = 0;
        let mut nodes_ignored: usize = 0;
        let mut peak_queue_size: usize = 0;
        let mut farthest_consensus: usize = 0;
        let mut last_constraint: u64 = 0;

        let offset_window = self.config.offset_window; // the window around the provided offset that we check
        let offset_compare_length = self.config.offset_compare_length; // the number of bases we check to figure out the real start

        let offsets: Vec<Option<usize>> = if self.config.auto_shift_offsets {
            // first figure out if we need to address the offsets for this re-run of consensus
            let mut min_offset = usize::MAX;
            let mut start_sequence_found = false;
            for offset in self.offsets.iter() {
                match offset {
                    Some(o) => min_offset = min_offset.min(*o),
                    None => start_sequence_found = true
                };
            }

            if !start_sequence_found {
                debug!("No start sequence detected, shifting all offsets by {min_offset}");
                self.offsets.iter()
                    .map(|o| {
                        let value = o.unwrap();
                        if value == min_offset {
                            // equal to minimum, flag as no offset
                            None
                        } else {
                            // past minimum, so down-shift it
                            Some(value - min_offset)
                        }
                    })
                    .collect::<Vec<Option<usize>>>()
            } else {
                self.offsets.clone()
            }
        } else {
            self.offsets.clone()
        };

        debug!("Offsets: {:?}", offsets);

        // build up the list of sizes where we need to activate one or more sequences
        let mut activate_points: HashMap<usize, Vec<usize>> = Default::default();
        let mut max_activate = 0;
        let mut initially_active: usize = 0;
        for (seq_index, offset) in offsets.iter().enumerate() {
            if let Some(last_offset) = offset {
                // we have an approx, calculate when it can get activated and add it to our list
                let activate_length = last_offset + offset_compare_length;
                let entry = activate_points.entry(activate_length).or_default();
                entry.push(seq_index);

                if activate_length > max_activate {
                    max_activate = activate_length;
                }
            } else {
                initially_active += 1;
            }
        }

        if initially_active == 0 {
            bail!("Must have at least one initial offset of None to see the consensus.");
        }
        
        let max_queue_size = self.config.max_queue_size;
        let max_capacity_per_size = self.config.max_capacity_per_size;
        let initial_size = self.sequences.iter().map(|s| s.len()).max().unwrap();
        let mut pqueue_tracker = PQueueTracker::with_capacity(initial_size, max_capacity_per_size);

        // settings and a vote buffer shared by `push()` and `summarize()`
        let mut step_context = StepContext {
            require_all: self.config.allow_early_termination,
            consensus_cost: self.config.consensus_cost,
            scratch: VoteScratch::default(),
        };

        let initial_node = ConsensusNode::new_root_node(
            &self.sequences, &offsets, &self.config, &mut step_context, terminal_scores,
        )?;

        // Max-heap ordered by queued priority: lower cost, then longer consensus, then earlier tie-breaker.
        let mut pqueue = PriorityQueue::new(self.consensus_cost());
        pqueue_tracker.insert(initial_node.consensus().len());
        pqueue.push(initial_node);

        // buffers that each search step reuses
        let mut extension_candidates: Vec<(u8, f64)> = vec![];
        let mut passing_candidates: Vec<u8> = vec![];

        // the way this will work is that we will eventually find one or more answers and anything worse will get drained off until no possibilities remain
        while !pqueue.is_empty() {
            // this just tracks how large our actual queue gets
            peak_queue_size = peak_queue_size.max(pqueue.len());

            // first, check if we need to restrict our pqueue threshold
            while (pqueue_tracker.len() > max_queue_size || last_constraint >= self.config.max_nodes_wo_constraint) && pqueue_tracker.threshold() < farthest_consensus {
                pqueue_tracker.increment_threshold();
                last_constraint = 0;
            }

            // get the top node and check if it's bad
            let top_node = pqueue.pop().unwrap();
            let top_priority = top_node.queued_priority();
            let top_len = top_node.consensus().len();
            pqueue_tracker.remove(top_len);

            trace!("Pop: {:?} => {:?}", top_priority, top_node.consensus);
        
            if top_priority.0.0 > results.maximum_error || top_priority.1 < pqueue_tracker.threshold() || pqueue_tracker.at_capacity(top_len) {
                trace!("\tIgnored");
                nodes_ignored += 1;
                continue;
            }

            // mark this one as explored and updated farthest if necessary
            farthest_consensus = farthest_consensus.max(top_len);
            nodes_explored += 1;
            last_constraint += 1;
            pqueue_tracker.process(top_len)?;

            // now check if this node has reached the end
            if top_node.summary.reached_end {
                consider_terminal(&top_node, &mut results)?;
            }

            // this fetches the list of options according to the WFA so far
            // NOTE: this CAN include the wildcard, but only if the wildcard is the only character
            // with 2+ candidates, the wildcard is removed, so something is left over
            top_node.candidates(self.config.wildcard, &mut extension_candidates);
            let max_observed = extension_candidates.iter().map(|&(_symbol, count)| count).max_by(|a, b| a.total_cmp(b))
                // if no observation, then just use min count
                .unwrap_or(self.config.min_count as f64);
            
            // the active threshold is the minimum of 1) the configured minimum count OR 2) the highest count we observed
            let active_threshold = (self.config.min_count as f64).min(max_observed);
            trace!("\tcandidates = {:?}", extension_candidates);

            // pull out the full list of viable candidates
            passing_candidates.clear();
            passing_candidates.extend(extension_candidates.iter()
                .filter(|&&(_symbol, count)| count >= active_threshold)
                .map(|&(symbol, _count)| symbol));

            if passing_candidates.is_empty() && top_len < max_activate {
                bail!("Encountered coverage gap: consensus is length {} with no candidates, but sequences activate at {}", top_len, max_activate);
            }
            // with no candidates, no extensions remain; this should just happen at the end

            // each child extends a clone; the last child takes `top_node` instead of a clone
            // one candidate extends in place without cloning
            // `top_node` is placed into an Option to pass the borrow checker
            let num_children = passing_candidates.len();
            let mut parent = Some(top_node);
            for (child_index, &symbol) in passing_candidates.iter().enumerate() {
                let mut new_node = if child_index + 1 == num_children {
                    parent.take().expect("parent is only taken on the last child")
                } else {
                    parent.as_ref().unwrap().clone()
                };
                new_node.push(&self.sequences, symbol, &mut step_context, &append_score_column)?;

                // check if we need to activate any strings
                let opt_activate_list = activate_points.get(&new_node.consensus().len());
                if let Some(activate_list) = opt_activate_list {
                    assert!(!activate_list.is_empty());
                    for &seq_index in activate_list.iter() {
                        new_node.activate_sequence(self.sequences[seq_index], seq_index, offset_window, offset_compare_length, &self.config)?;
                    }
                    // summarize once after every read at this length has been activated
                    new_node.summarize(&self.sequences, &mut step_context);
                }

                // get the new cost and put it in the queue
                trace!("\tPush {} => {:?}", pqueue.next_id(), new_node.consensus);
                pqueue_tracker.insert(new_node.consensus().len());
                pqueue.push(new_node);
            }
        }

        assert_eq!(pqueue_tracker.len(), 0);

        // sort these by the sequence so we always have a fixed order
        results.retained.sort_by(|c1, c2| c1.sequence().cmp(c2.sequence()));

        debug!("nodes_explored: {nodes_explored}");
        debug!("nodes_ignored: {nodes_ignored}");
        debug!("peak_queue_size: {peak_queue_size}");
        Ok(results.retained)
    }

    // getters
    pub fn sequences(&self) -> &[&'a [u8]] {
        &self.sequences
    }

    pub fn consensus_cost(&self) -> ConsensusCost {
        self.config.consensus_cost
    }

    pub fn alphabet(&self) -> &HashSet<u8> {
        &self.alphabet
    }
}

/// The best terminal cost also bounds exploration using the DWFA prefix cost.
struct TerminalResults {
    maximum_error: usize,
    retained: Vec<Consensus>,
}

impl TerminalResults {
    fn retain(
        &mut self,
        score: usize,
        max_return_size: usize,
        build: impl FnOnce() -> Result<Consensus, SimpleError>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if score < self.maximum_error {
            self.maximum_error = score;
            self.retained.clear();
        }
        if score <= self.maximum_error && self.retained.len() < max_return_size {
            self.retained.push(build()?);
        }
        Ok(())
    }
}

/// Cached state for the next search step.
/// `push()` refreshes it for each extension.
/// `new_root_node()` refreshes it through `summarize()`.
/// After one or more activations, the caller refreshes it once with `summarize()`.
/// A node does not change once it enters the search heap, so this cache remains valid.
#[derive(Clone, Debug, Default)]
struct StepSummary {
    /// The extension votes of all active DWFAs, each with weight 1.0.
    /// Symbols are ordered by their first vote as the reads are traversed.
    votes: ExtensionVotes,
    /// Whether this node may be finalized. When `StepContext::require_all` is false, any active DWFA at its
    /// baseline end is sufficient; when true, every DWFA must be at its baseline end. An inactive DWFA is not
    /// at the end, while a capped DWFA contributes the `require_all` value.
    reached_end: bool,
    /// Sum of the active DWFA costs under `StepContext::consensus_cost`.
    total_cost: usize,
}

/// The settings of `StepSummary` and a scratch buffer, for one consensus run.
struct StepContext {
    /// Require every read to reach its end before finalization, as configured by `allow_early_termination`.
    require_all: bool,
    /// Cost model used for `StepSummary::total_cost` and priority-queue ordering.
    consensus_cost: ConsensusCost,
    /// Reusable workspace for collecting a DWFA's extension vote.
    scratch: VoteScratch,
}

impl StepContext {
    /// The cost of one DWFA with this edit distance
    /// # Arguments
    /// * `edit_distance` - the edit distance of one DWFA
    #[inline(always)]
    fn cost(&self, edit_distance: usize) -> usize {
        match self.consensus_cost {
            ConsensusCost::L1Distance => edit_distance,
            ConsensusCost::L2Distance => edit_distance.pow(2)
        }
    }

    /// Whether one DWFA counts as at the end for `StepSummary::reached_end`.
    /// A capped read must not finish the consensus early, and must not block when every read has to finish,
    /// so it gives the `require_all` value.
    #[inline(always)]
    fn at_end(&self, dwfa: &DWFALite, baseline: &[u8]) -> bool {
        if dwfa.state() == DWFALiteState::ExceededEditDistanceLimit {
            self.require_all
        } else {
            dwfa.reached_baseline_end(baseline)
        }
    }

    /// `StepSummary::reached_end` from the any and all values over the reads
    #[inline(always)]
    fn reached_end(&self, any_at_end: bool, all_at_end: bool) -> bool {
        if self.require_all { all_at_end } else { any_at_end }
    }
}

/// Wrapper for a node containing a partial consensus as well as the DWFA tracking for that node.
#[derive(Clone, Debug)]
struct ConsensusNode<S> {
    /// The consensus sequence so far
    consensus: Vec<u8>,
    /// The DWFAs that are tracked for each sequence; these are only None if they have not been activated yet due to an offset
    dwfas: Vec<Option<DWFALite>>,
    /// Rank in the search heap: lower cost, then longer consensus, then earlier tie-breaker. Set by `cache_priority`.
    queued_priority: NodePriority,
    /// What the next search step needs from `dwfas`. Set by `push()`, `new_root_node()`, and one `summarize()` after activation.
    summary: StepSummary,
    /// Incremental terminal scores for this branch, or `()` for native WFA finalization.
    terminal_scores: S,
}

// The following Eq and Ord implementations are specifically for the Binary Heap.
// They should not be used for actual node comparison.

impl<S> PartialEq for ConsensusNode<S> {
    fn eq(&self, other: &Self) -> bool {
        self.queued_priority == other.queued_priority
    }
}

// `PartialEq` only compares the priority, which is a total order
impl<S> Eq for ConsensusNode<S> {}

impl<S> PartialOrd for ConsensusNode<S> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<S> Ord for ConsensusNode<S> {
    fn cmp(&self, other: &Self) -> Ordering {
        // BinaryHeap pops the greatest node. queued_priority already ranks a lower edit distance,
        // then a longer consensus, then an earlier tie-breaker.
        self.queued_priority.cmp(&other.queued_priority)
    }
}

impl<S> ConsensusNode<S> {
    /// Constructor for a new consensus search root node
    /// # Arguments
    /// * `sequences` - the baseline reads, parallel to `offsets`
    /// * `offsets` - a set of offsets into the sequences where the approximate starts are
    /// * `config` - consensus configuration
    /// * `context` - summary settings; the returned node is already summarized
    /// * `terminal_scores` - terminal score state for the empty consensus
    /// # Errors
    /// * if DWFA construction fails
    /// * if `sequences` and `offsets` have different lengths
    fn new_root_node(
        sequences: &[&[u8]],
        offsets: &[Option<usize>],
        config: &CdwfaConfig,
        context: &mut StepContext,
        terminal_scores: S,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if sequences.len() != offsets.len() {
            bail!("Sequence and offset counts must match.");
        }
        let mut dwfas = Vec::with_capacity(offsets.len());
        for (sequence, offset) in sequences.iter().zip(offsets.iter()) {
            let dwfa = if offset.is_some() {
                // we have an offset, so do not create a map
                None
            } else {
                // active from the start; cap is fixed from this read's length
                Some(DWFALite::new(config.dwfa_lite_config_for(sequence.len())?))
            };
            dwfas.push(dwfa);
        }

        let active_count = dwfas.iter().filter(|d| d.is_some()).count();
        if active_count == 0 {
            bail!("Root consensus node has no initially active sequences.");
        }

        let mut node = Self {
            consensus: vec![],
            dwfas,
            queued_priority: (Reverse(0), 0, Reverse(0)),
            summary: StepSummary::default(),
            terminal_scores,
        };
        node.summarize(sequences, context);
        Ok(node)
    }

    /// This will activate the DWFA for a particular sequence.
    /// 
    /// NOTE: This does *not* summarize the node.
    /// Users must separately call `summarize()` after activating all sequences at this consensus length.
    /// # Arguments
    /// * `sequence` - the sequence that is getting activated
    /// * `seq_index` - the sequence index, which will map to the DWFAs
    /// * `offset_window` - the window we are searching for a best match
    /// * `offset_compare_length` - the amount of bases we are comparing
    /// * `config` - consensus configuration, including wildcard, early termination, and the edit-distance cap
    /// # Errors
    /// * if DWFA construction or the initial update fails
    fn activate_sequence(&mut self, sequence: &[u8], seq_index: usize, offset_window: usize, offset_compare_length: usize, config: &CdwfaConfig) -> Result<(), Box<dyn std::error::Error>> {
        // make sure everything is currently inactive
        assert!(self.dwfas[seq_index].is_none());

        // build up the search space
        let con_len = self.consensus.len();

        // we need to truncate the compare length if it is longer than the sequence
        let offset_compare_length = offset_compare_length.min(sequence.len());
        
        // figure out how far back we search
        let start_delta = offset_window + offset_compare_length; // we search from current back to the offset_window
        let start_position = con_len.saturating_sub(start_delta);
        let end_position = con_len.saturating_sub(offset_compare_length);

        // figure out which offset has the best score; assume the middle of the offset window is the best
        let mut best_offset = con_len.saturating_sub(offset_compare_length + offset_window / 2);
        let mut min_ed = crate::sequence_alignment::wfa_ed_config(&self.consensus[best_offset..], &sequence[0..offset_compare_length], false, config.wildcard);
        
        // now check all the rest around this position
        for p in start_position..end_position {
            let ed = crate::sequence_alignment::wfa_ed_config(&self.consensus[p..], &sequence[0..offset_compare_length], false, config.wildcard);
            if ed < min_ed {
                min_ed = ed;
                best_offset = p;
            }
        }

        // now set up the DWFA with the best offset and a cap from this read's length
        let mut new_dwfa = DWFALite::new(config.dwfa_lite_config_for(sequence.len())?);
        new_dwfa.set_offset(best_offset);
        new_dwfa.update(sequence, &self.consensus)?;
        self.dwfas[seq_index] = Some(new_dwfa);

        Ok(())
    }

    /// Extends the consensus, active wavefronts, and terminal scores with the same symbol.
    /// A single next symbol returned by `update()` votes directly; otherwise the DWFA calls `vote()`.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `symbol` - the symbol to extend each DWFA with
    /// * `context` - summary settings and vote scratch space
    /// * `append_score_column` - advances the terminal scores for this branch
    /// # Errors
    /// * if any of the extensions fail
    fn push(
        &mut self,
        baseline_sequences: &[&[u8]],
        symbol: u8,
        context: &mut StepContext,
        append_score_column: &impl Fn(&mut S, u8),
    ) -> Result<(), Box<dyn std::error::Error>> {
        // update the consensus
        self.consensus.push(symbol);

        // clear the votes in preparation for the new symbols
        let votes = &mut self.summary.votes;
        votes.clear();

        // trackers going into the summary
        let mut total_cost = 0;
        let mut any_at_end = false;
        let mut all_at_end = true;
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter_mut()) {
            let at_end = match opt_dwfa.as_mut() {
                Some(dwfa) => {
                    let (edit_distance, next_symbol) = dwfa.update(baseline, &self.consensus)?;
                    match next_symbol {
                        Some(symbol) => votes.add_single(symbol, 1.0),
                        None => dwfa.vote(baseline, &self.consensus, 1.0, votes, &mut context.scratch),
                    }
                    total_cost += context.cost(edit_distance);
                    context.at_end(dwfa, baseline)
                },
                None => false,
            };
            any_at_end |= at_end;
            all_at_end &= at_end;
        }
        self.summary.total_cost = total_cost;
        self.summary.reached_end = context.reached_end(any_at_end, all_at_end);
        append_score_column(&mut self.terminal_scores, symbol);
        Ok(())
    }

    /// Computes votes, total cost, and `reached_end` from the current DWFAs.
    /// Called by `new_root_node()`, and should be called after activating new sequences.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `context` - the summary settings
    fn summarize(&mut self, baseline_sequences: &[&[u8]], context: &mut StepContext) {
        let votes = &mut self.summary.votes;
        votes.clear();
        let mut total_cost = 0;
        let mut any_at_end = false;
        let mut all_at_end = true;
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter()) {
            let at_end = match opt_dwfa.as_ref() {
                Some(dwfa) => {
                    dwfa.vote(baseline, &self.consensus, 1.0, votes, &mut context.scratch);
                    total_cost += context.cost(dwfa.edit_distance());
                    context.at_end(dwfa, baseline)
                },
                None => false,
            };
            any_at_end |= at_end;
            all_at_end &= at_end;
        }
        self.summary.total_cost = total_cost;
        self.summary.reached_end = context.reached_end(any_at_end, all_at_end);
    }

    /// Finalizes all of the contained DWFAs
    /// # Errors
    /// * if any of the individual finalizes fail
    fn finalize(&mut self, baseline_sequences: &[&[u8]]) -> Result<(), Box<dyn std::error::Error>> {
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter_mut()) {
            if let Some(dwfa) = opt_dwfa.as_mut() {
                #[cfg(test)]
                NATIVE_FINALIZATIONS.with(|count| count.set(count.get() + 1));
                dwfa.finalize(baseline, &self.consensus)?;
            } else {
                bail!("Finalize called on DWFA that was never initialized.");
            }
        }
        Ok(())
    }

    /// Returns a vec of all of the scores for this node
    fn costs(&self, consensus_cost: ConsensusCost) -> Vec<usize> {
        self.dwfas.iter()
            .map(|opt_d| {
                if let Some(d) = opt_d.as_ref() {
                    match consensus_cost {
                        ConsensusCost::L1Distance => d.edit_distance(),
                        ConsensusCost::L2Distance => d.edit_distance().pow(2)
                    }
                } else {
                    // it's not initialized, here we will just report 0
                    0
                }
            })
            .collect()
    }

    /// Labels each read. `None` when every read was used.
    /// `Some` when any read was not Included.
    fn assignments(&self) -> Option<Vec<ConsensusAssignment>> {
        let labels: Vec<ConsensusAssignment> = self.dwfas.iter()
            .map(|opt_dwfa| {
                if opt_dwfa.as_ref().is_some_and(|dwfa| dwfa.state() == DWFALiteState::ExceededEditDistanceLimit) {
                    ConsensusAssignment::EditDistanceLimit
                } else {
                    ConsensusAssignment::Included
                }
            })
            .collect();

        // if ALL labels are included, then we can safely send back None
        if labels.iter().all(|label| *label == ConsensusAssignment::Included) {
            None
        } else {
            // otherwise, we need to send back the labels
            Some(labels)
        }
    }

    /// Returns the total score for the node
    fn total_cost(&self, consensus_cost: ConsensusCost) -> usize {
        // same as summing `costs()`, without the allocation; a DWFA that is not initialized costs 0
        self.dwfas.iter()
            .flatten()
            .map(|d| match consensus_cost {
                ConsensusCost::L1Distance => d.edit_distance(),
                ConsensusCost::L2Distance => d.edit_distance().pow(2)
            })
            .sum()
    }

    /// Writes this node's extension candidates into `out`.
    /// When there are two or more candidates, `wildcard` is removed. One candidate is kept as it is,
    /// including when that candidate is the wildcard.
    /// # Arguments
    /// * `wildcard` - an optional wildcard symbol
    /// * `out` - output buffer; it is cleared first
    fn candidates(&self, wildcard: Option<u8>, out: &mut Vec<(u8, f64)>) {
        self.summary.votes.candidates(wildcard, out);
    }

    // getters
    fn consensus(&self) -> &[u8] {
        &self.consensus
    }

    /// Rank cached for the search heap: lower cost, then longer consensus, then earlier tie-breaker.
    fn queued_priority(&self) -> NodePriority {
        self.queued_priority
    }
}

impl<S> CachePriority for ConsensusNode<S> {
    /// Caches this node's search priority.
    /// Rank is 1) lowest cost, 2) consensus length, and 3) earlier tie-breaker.
    fn cache_priority(&mut self, consensus_cost: ConsensusCost, tie_breaker: u64) {
        // the summary has the cost from the config, which the queue also uses
        debug_assert_eq!(self.summary.total_cost, self.total_cost(consensus_cost));
        self.queued_priority = (
            Reverse(self.summary.total_cost),
            self.consensus.len(),
            Reverse(tie_breaker),
        );
    }
}

#[cfg(test)]
std::thread_local! {
    static NATIVE_FINALIZATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BinaryHeap;

    use crate::cdwfa_config::CdwfaConfigBuilder;

    fn with_native_finalization_count<T>(run: impl FnOnce() -> T) -> (T, usize) {
        NATIVE_FINALIZATIONS.with(|count| {
            count.set(0);
            (run(), count.get())
        })
    }

    #[test]
    fn test_global_scores_match_length_divergent_reads() {
        let reads: [&[u8]; 4] = [b"ACGT", b"AC", b"A", b"ACG"];
        let mut dwfa = ConsensusDWFA::default();
        for read in reads {
            dwfa.add_sequence(read).unwrap();
        }

        let (legacy, legacy_finalizations) =
            with_native_finalization_count(|| dwfa.consensus().unwrap());
        assert_eq!(legacy.len(), 2);
        assert_eq!(legacy[0].sequence(), b"AC");
        assert_eq!(legacy[1].sequence(), b"ACG");
        // The length-one, length-two, and length-three terminals each finalize four reads.
        assert_eq!(legacy_finalizations, 12);

        let (global, global_finalizations) =
            with_native_finalization_count(|| dwfa.consensus_with_global_scores().unwrap());
        assert_eq!(global, legacy);
        assert_eq!(global_finalizations, 0);
    }

    #[test]
    fn test_global_scores_match_an_empty_read_and_terminal_root() {
        let mut dwfa = ConsensusDWFA::default();
        dwfa.add_sequence(b"ACGT").unwrap();
        dwfa.add_sequence(b"").unwrap();

        let legacy = dwfa.consensus().unwrap();
        assert_eq!(legacy[0].sequence(), b"");
        assert_eq!(legacy[0].scores(), &[4, 0]);
        assert_eq!(dwfa.consensus_with_global_scores().unwrap(), legacy);
    }

    #[test]
    fn test_global_scores_preserve_tied_siblings_and_the_return_cap() {
        let mut first = [b'A'; 65];
        first[64] = 0xff;
        let mut second = first;
        second[64] = 0;
        for max_return_size in [1, 2] {
            let config = CdwfaConfigBuilder::default()
                .max_return_size(max_return_size)
                .build()
                .unwrap();
            let mut dwfa = ConsensusDWFA::with_config(config).unwrap();
            dwfa.add_sequence(&first).unwrap();
            dwfa.add_sequence(&second).unwrap();

            let legacy = dwfa.consensus().unwrap();
            assert_eq!(legacy.len(), max_return_size);
            // The first queued branch wins a capped tie even though it sorts last.
            assert_eq!(legacy.last().unwrap().sequence(), &first[..]);
            assert_eq!(dwfa.consensus_with_global_scores().unwrap(), legacy);
        }
    }

    #[test]
    fn test_global_scores_fall_back_for_unsupported_configurations() {
        let configs = [
            CdwfaConfig { consensus_cost: ConsensusCost::L2Distance, ..Default::default() },
            CdwfaConfig { wildcard: Some(b'T'), ..Default::default() },
            CdwfaConfig { allow_early_termination: true, ..Default::default() },
            CdwfaConfig { max_edit_distance: Some(3), ..Default::default() },
            CdwfaConfig { max_edit_distance_fraction: Some(0.25), ..Default::default() },
            CdwfaConfig { min_edit_distance: Some(3), ..Default::default() },
        ];
        for config in configs {
            let mut dwfa = ConsensusDWFA::with_config(config.clone()).unwrap();
            for read in [b"ACGTACGTACGT", b"ACGTACGTACGT", b"ACGTACGTACGT", b"TTTTTTTTTTTT"] {
                dwfa.add_sequence(read).unwrap();
            }

            let (legacy, legacy_finalizations) =
                with_native_finalization_count(|| dwfa.consensus().unwrap());
            let (global, global_finalizations) =
                with_native_finalization_count(|| dwfa.consensus_with_global_scores().unwrap());
            assert_eq!(global, legacy, "{config:?}");
            assert!(legacy_finalizations > 0, "{config:?}");
            assert_eq!(global_finalizations, legacy_finalizations, "{config:?}");
        }
    }

    /// Checks that a lower cost wins, and equal costs prefer an earlier tie-breaker.
    #[test]
    fn test_heap_orders_by_priority_then_id() {
        let sequences: [&[u8]; 1] = [b""];
        let offsets = vec![None];
        let config = CdwfaConfig::default();
        let mut context = StepContext {
            require_all: config.allow_early_termination,
            consensus_cost: config.consensus_cost,
            scratch: VoteScratch::default(),
        };
        let mut earlier = ConsensusNode::new_root_node(&sequences, &offsets, &config, &mut context, ()).unwrap();
        let mut later = ConsensusNode::new_root_node(&sequences, &offsets, &config, &mut context, ()).unwrap();
        let mut cheaper = ConsensusNode::new_root_node(&sequences, &offsets, &config, &mut context, ()).unwrap();
        earlier.queued_priority = (Reverse(1), 0, Reverse(0));
        later.queued_priority = (Reverse(1), 0, Reverse(1));
        cheaper.queued_priority = (Reverse(0), 0, Reverse(2));

        let mut heap = BinaryHeap::new();
        heap.push(later);
        heap.push(earlier);
        heap.push(cheaper);

        assert_eq!(heap.pop().unwrap().queued_priority().2.0, 2);
        assert_eq!(heap.pop().unwrap().queued_priority().2.0, 0);
        assert_eq!(heap.pop().unwrap().queued_priority().2.0, 1);
    }

    /// Checks that one exact input is returned as the consensus with score 0.
    #[test]
    fn test_single_sequence() {
        let sequence = b"ACGTACGTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();
        consensus_dwfa.add_sequence(sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![Consensus {
            sequence: sequence.to_vec(),
            consensus_cost: ConsensusCost::L1Distance,
            scores: vec![0],
            assignments: None,
        }]);
    }

    /// Checks that two sequences one edit apart both win, ordered lexicographically.
    #[test]
    fn test_dual_sequence() {
        let sequence = b"ACGTACGTACGT";
        let sequence2 = b"ACGTACCTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence2).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is a tie between the two we put it, each 1 ED total away
        let consensus = consensus_dwfa.consensus().unwrap();

        // sequence2 is alphabetically before sequence 1, so it will come first
        assert_eq!(consensus, vec![
            Consensus {
                sequence: sequence2.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![1, 0],
                assignments: None,
            },
            Consensus {
                sequence: sequence.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![0, 1],
                assignments: None,
            },
        ]);
    }

    /// Checks that the majority sequence is the unique consensus when one read differs by one edit.
    #[test]
    fn test_trio_sequence() {
        let sequence = b"ACGTACGTACGT";
        let sequence2 = b"ACGTACCTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();

        // add the first sequence twice, it will be the consensus
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence2).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found sequence
        let consensus = consensus_dwfa.consensus().unwrap();

        // sequence2 is alphabetically before sequence 1, so it will come first
        assert_eq!(consensus, vec![
            Consensus {
                sequence: sequence.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![0, 0, 1],
                assignments: None,
            }
        ]);
    }

    /// Checks that three indel-divergent reads recover a single intermediate consensus.
    #[test]
    fn test_complicated() {
        let expected_consensus = b"ACGTACGTACGT";
        let sequences = [
            // b"ACGTACG-TACGT"
            // b"AC-TACGGTACGT"
            b"ACTACGGTACGT",
            // b"ACGTAAG-TCCGT"
            b"ACGTAAGTCCGT",
            // b"AAGTACG-TACGT"
            b"AAGTACGTACGT"
        ];
        
        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::default();
        for &sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
    }

    /// Checks that wildcard bases are ignored for the alphabet and for edit scoring.
    #[test]
    fn test_wildcards() {
        let expected_consensus = b"ACGTACGTACGT";
        let sequences = [
            // C insertion and tail missing
            b"ACGTACCGT****".to_vec(),
            // C>T and head+tail missing
            b"**GTATGTAC**".to_vec(),
            // just head missing
            b"****ACGTACGT".to_vec()
        ];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
        assert_eq!(&[1, 1, 0], consensus[0].scores());
    }

    /// Checks that a lone wildcard remains a candidate when every read nominates it.
    #[test]
    fn test_all_wildcards() {
        let _expected_consensus = b"ACGTACGTACGT";
        let actual_consensus = b"*CGTACG*ACG*";
        let sequences = [
            // 1 ins
            b"*CGTAACG*ACG*".to_vec(),
            // nothing other than *
            b"*CGTACG*ACG*".to_vec(),
            // C>T
            b"*CGTACG*ATG*".to_vec()
        ];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(actual_consensus, consensus[0].sequence());
        assert_eq!(&[1, 0, 1], consensus[0].scores());
    }

    /// Checks that early termination recovers a full consensus from prefix-truncated reads.
    #[test]
    fn test_allow_early_termination_costs() {
        let expected_consensus = b"ACGT";

        // first verify baseline is not generating the full thing due to error
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for i in 1..(expected_consensus.len()+1) {
            consensus_dwfa.add_sequence(&expected_consensus[0..i]).unwrap();
        }
        // this first approach generated multiple possible ones that are in the middle
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, [
            Consensus { sequence: vec![65, 67], consensus_cost: ConsensusCost::L1Distance, scores: vec![1, 0, 1, 2], assignments: None },
            Consensus { sequence: vec![65, 67, 71], consensus_cost: ConsensusCost::L1Distance, scores: vec![2, 1, 0, 1], assignments: None }
        ]);

        // second, verify that allowing early termination fixes it to the original result
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .allow_early_termination(true)
                .build().unwrap()
        ).unwrap();
        for i in 1..(expected_consensus.len()+1) {
            consensus_dwfa.add_sequence(&expected_consensus[0..i]).unwrap();
        }
        // this first approach generated multiple possible ones that are in the middle
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, [
            Consensus { sequence: expected_consensus.to_vec(), consensus_cost: ConsensusCost::L1Distance, scores: vec![0; 4], assignments: None },
        ]);
    }

    /// Checks that offset windows place delayed starts so every read scores 0.
    #[test]
    fn test_offset_windows() {
        let expected_consensus = b"ACGTACGTACGTACGT";
        let sequences = [
            // no offset
            b"ACGTACGTACGTACGT".to_vec(),
            // offset 4
            b"ACGTACGTACGT".to_vec(),
            // offset 6
            b"GTACGTACGT".to_vec(),
        ];

        // make one offset exactly match and another be right on our boundary (1-shift here)
        let offsets = [None, Some(4), Some(7)];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .offset_window(1) // both are exactly one off
                .offset_compare_length(4) // semi-arbitrary, just needs to be short here given our seq-len
                .build().unwrap()
        ).unwrap();
        for (sequence, offset) in sequences.iter().zip(offsets.iter().cloned()) {
            consensus_dwfa.add_sequence_offset(sequence, offset).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
        assert_eq!(&[0, 0, 0], consensus[0].scores());
        assert_eq!(consensus_dwfa.consensus_with_global_scores().unwrap(), consensus);
    }

    /// Checks that a start far past activation leaves a DWFA uninitialized and fails finalize.
    #[test]
    fn test_offset_gap_err() {
        let sequences = [
            // no offset
            b"ACGTACGTACGTACGT".to_vec(),
            // massive offset
            b"ACGTACGTACGTACGT".to_vec(),
        ];

        // one starts and one is way out there
        let offsets = [None, Some(1000)];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .offset_window(1) // both are exactly one off
                .offset_compare_length(4) // semi-arbitrary, just needs to be short here given our seq-len
                .build().unwrap()
        ).unwrap();
        for (sequence, offset) in sequences.iter().zip(offsets.iter().cloned()) {
            consensus_dwfa.add_sequence_offset(sequence, offset).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus_err = consensus_dwfa.consensus();
        assert!(consensus_err.is_err());
        let error = consensus_err.err().unwrap();
        assert_eq!(error.to_string(), "Finalize called on DWFA that was never initialized.");
        let global_error = consensus_dwfa.consensus_with_global_scores().unwrap_err();
        assert_eq!(global_error.to_string(), error.to_string());
    }

    /// Checks that a hard edit-distance cap freezes an outlier and labels it in the result.
    #[test]
    fn test_edit_distance_cap_ignores_outlier() {
        let good = b"ACGTACGTACGT";
        let bad = b"TTTTTTTTTTTT";
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .max_edit_distance(Some(3))
                .build().unwrap()
        ).unwrap();
        consensus_dwfa.add_sequence(good).unwrap();
        consensus_dwfa.add_sequence(good).unwrap();
        consensus_dwfa.add_sequence(good).unwrap();
        consensus_dwfa.add_sequence(bad).unwrap();

        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(consensus[0].sequence(), good);
        assert_eq!(consensus[0].scores(), &[0, 0, 0, 3]);
        assert_eq!(consensus[0].assignments(), Some(&[
            ConsensusAssignment::Included,
            ConsensusAssignment::Included,
            ConsensusAssignment::Included,
            ConsensusAssignment::EditDistanceLimit,
        ][..]));
    }
}