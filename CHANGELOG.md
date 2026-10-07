# Unreleased
## Performance
* Integrate upstream's binary-heap consensus queues and dual-consensus benchmarks.
* Remove the superseded queue hashing and cached consensus-prefix fingerprints.
* Retain reusable candidate storage, allocation-free total-cost calculation, in-place wavefront growth, and maximum-reach tracking during extension.

## Features
* Add optional absolute and fractional per-read edit-distance caps with an optional minimum cap. Reads exceeding a cap stop voting and expose an explicit exclusion assignment; caps remain disabled by default.

## API changes
* `DWFALite::new` takes `DWFALiteConfig`, with lifecycle available through `DWFALiteState`.
* `Consensus::new` takes an optional assignment vector and returns a `Result`.
* Dual-consensus `assignments()` replaces `is_consensus1()` with explicit allele, equal-score, and edit-distance-limit states.

## Bugs
* Preserve exact global endpoint scoring and finalization sealing while integrating edit-distance limits, including caps reached during finalization.

# v0.4.4
## Bugs
* Fixed a panic caused by input sequences shorter than the offset compare length

# v0.4.3
## Bugs
* Fixed a panic caused by an unchecked `unwrap()` in a `trace!` statement

# v0.4.2
Initial release
