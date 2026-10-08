# Unreleased
## Performance
* Integrate upstream's binary-heap consensus queues and dual-consensus benchmarks.
* Remove the superseded queue hashing and cached consensus-prefix fingerprints.
* Retain reusable candidate storage, allocation-free total-cost calculation, in-place wavefront growth, and maximum-reach tracking during extension.
* Enqueue consensus children directly and move the parent into the last child, avoiding a temporary child vector and one deep clone per branch.
* Cache exact baseline and other-coordinate maxima while preserving DWFA equality, hashing, offsets, caps, and finalization semantics.
* Extend generated-read benchmarks with 10 kb single-consensus reads and a branching fixture. On Apple arm64 with Rust 1.85.1/default-feature release LTO, the retained ownership-plus-maxima variant improved the 1 kb/30-read/1% and branching cases by 10.0% and 16.7%, respectively (median of three paired ratios), with no regression above 3% across 15 single and 12 dual cases. Rightmost-interval scanning was evaluated but not retained because seven full-matrix rows exceeded that guard; CNBP was excluded.

## Features
* Add optional absolute and fractional per-read edit-distance caps with an optional minimum cap. Reads exceeding a cap stop voting and expose an explicit exclusion assignment; caps remain disabled by default.
* Add a default-enabled `logging` Cargo feature. Disabling default features removes logging calls, diagnostic search statistics, and the normal `log` dependency without removing consensus functionality.

## API changes
* `DWFALite::new` takes `DWFALiteConfig`, with lifecycle available through `DWFALiteState`.
* `Consensus::new` takes an optional assignment vector and returns a `Result`.
* Dual-consensus `assignments()` replaces `is_consensus1()` with explicit allele, equal-score, and edit-distance-limit states.

## Bugs
* Preserve exact global endpoint scoring and finalization sealing while integrating edit-distance limits, including caps reached during finalization.
* Make dual-consensus trace formatting byte-safe so enabling logging cannot panic or return UTF-8 conversion errors for non-UTF-8 sequences.

## Packaging
* Add package description, repository, documentation, license-file metadata, and a declared Rust 1.85.1 minimum.
* Restrict published contents to source, tests and fixtures, benchmarks, examples, and package documentation; exclude internal research artifacts and CNBP datasets.
* Exercise logging-enabled and logging-disabled configurations on stable Rust and the minimum supported version in CI.

# v0.4.4
## Bugs
* Fixed a panic caused by input sequences shorter than the offset compare length

# v0.4.3
## Bugs
* Fixed a panic caused by an unchecked `unwrap()` in a `trace!` statement

# v0.4.2
Initial release
