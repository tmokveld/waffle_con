use waffle_con::dual_consensus::{DualConsensus, DualConsensusDWFA, SequenceAssignment};

#[cfg(feature = "logging")]
struct SinkLogger;

#[cfg(feature = "logging")]
impl log::Log for SinkLogger {
    fn enabled(&self, _: &log::Metadata<'_>) -> bool {
        true
    }

    fn log(&self, record: &log::Record<'_>) {
        use std::io::Write;
        std::io::sink().write_fmt(*record.args()).unwrap();
    }

    fn flush(&self) {}
}

fn non_utf8_consensuses() -> [Vec<DualConsensus>; 2] {
    let sequence = [0xff];
    let mut single = DualConsensusDWFA::default();
    for _ in 0..3 {
        single.add_sequence(&sequence).unwrap();
    }
    let single_results = single.consensus().unwrap();
    assert_eq!(single_results.len(), 1);
    assert_eq!(single_results[0].consensus1().sequence(), &[0xff]);
    assert_eq!(single_results[0].consensus1().scores(), &[0, 0, 0]);
    assert!(!single_results[0].is_dual());

    let allele1 = [0xfe, b'A'];
    let allele2 = [0xff, b'A'];
    let mut dual = DualConsensusDWFA::default();
    for allele in [&allele1, &allele2] {
        for _ in 0..3 {
            dual.add_sequence(allele).unwrap();
        }
    }
    let dual_results = dual.consensus().unwrap();
    assert_eq!(dual_results.len(), 1);
    let result = &dual_results[0];
    assert_eq!(result.consensus1().sequence(), &[0xfe, b'A']);
    assert_eq!(result.consensus2().unwrap().sequence(), &[0xff, b'A']);
    assert_eq!(result.consensus1().scores(), &[0, 0, 0]);
    assert_eq!(result.consensus2().unwrap().scores(), &[0, 0, 0]);
    assert_eq!(result.assignments(), &[
        SequenceAssignment::Consensus1,
        SequenceAssignment::Consensus1,
        SequenceAssignment::Consensus1,
        SequenceAssignment::Consensus2,
        SequenceAssignment::Consensus2,
        SequenceAssignment::Consensus2,
    ]);
    [single_results, dual_results]
}

#[test]
fn logging_preserves_non_utf8_consensus_results() {
    // This integration test runs in its own process, isolating the global logger.
    #[cfg(feature = "logging")]
    {
        static LOGGER: SinkLogger = SinkLogger;
        log::set_logger(&LOGGER).unwrap();
        log::set_max_level(log::LevelFilter::Off);
    }

    let without_logging = non_utf8_consensuses();
    #[cfg(feature = "logging")]
    {
        log::set_max_level(log::LevelFilter::Trace);
        let with_logging = non_utf8_consensuses();
        assert_eq!(with_logging, without_logging);
        // DualConsensus equality does not include the per-read score arrays.
        for (logged, silent) in with_logging.iter().zip(&without_logging) {
            assert_eq!(logged[0].scores1(), silent[0].scores1());
            assert_eq!(logged[0].scores2(), silent[0].scores2());
        }
        log::set_max_level(log::LevelFilter::Off);
    }
    #[cfg(not(feature = "logging"))]
    let _ = without_logging;
}
