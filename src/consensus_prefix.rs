use rustc_hash::FxHasher;
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::Deref;

/// Append-only bytes with a constant-size hash; equality always compares bytes.
#[derive(Clone, Default)]
pub(crate) struct ConsensusPrefix {
    sequence: Vec<u8>,
    fingerprint: u64,
}

impl ConsensusPrefix {
    pub(crate) fn push(&mut self, symbol: u8) {
        self.sequence.push(symbol);
        let mut hasher = FxHasher::default();
        hasher.write_u64(self.fingerprint);
        hasher.write_u8(symbol);
        self.fingerprint = hasher.finish();
    }
}

impl Deref for ConsensusPrefix {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.sequence
    }
}

impl fmt::Debug for ConsensusPrefix {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.sequence.fmt(formatter)
    }
}

impl PartialEq for ConsensusPrefix {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence
    }
}

impl Eq for ConsensusPrefix {}

impl PartialOrd for ConsensusPrefix {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ConsensusPrefix {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sequence.cmp(&other.sequence)
    }
}

impl Hash for ConsensusPrefix {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence.len().hash(state);
        self.fingerprint.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hash::FxHashMap;

    #[test]
    fn test_prefix_hash_tracks_forks() {
        let mut first = ConsensusPrefix::default();
        first.push(0);
        first.push(255);
        let mut second = first.clone();
        first.push(0);
        second.push(255);
        assert_eq!(first.cmp(&second), first[..].cmp(&second[..]));

        let mut entries = FxHashMap::default();
        entries.insert(first.clone(), 1);
        entries.insert(second.clone(), 2);
        let mut rebuilt = ConsensusPrefix::default();
        for &symbol in first.iter() {
            rebuilt.push(symbol);
        }
        assert_eq!(entries.get(&rebuilt), Some(&1));
        assert_eq!(entries.insert(rebuilt, 3), Some(1));
        assert_eq!(entries.len(), 2);
        assert_eq!(entries.get(&first), Some(&3));
        assert_eq!(entries.get(&second), Some(&2));
    }
}
