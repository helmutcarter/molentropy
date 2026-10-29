//! Estimates are bit-for-bit reproducible: the same input gives the same
//! floating-point result on every run and with any number of threads.

mod common;

use std::f64::consts::TAU;

use common::Rng;
use molentropy::{
    CoordinateMetric, EntropyOptions, calculate_entropy, estimate_coordinate_entropy,
    estimate_coordinate_network_entropy, estimate_network_entropy,
};

const FRAMES: usize = 200;

/// Runs `estimate` on a fresh pool of `threads` threads.
fn with_threads<T: Send>(threads: usize, estimate: impl FnOnce() -> T + Send) -> T {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(estimate)
}

/// Alternating linear and periodic coordinates.
fn system(dimensions: usize) -> (Vec<Vec<f64>>, Vec<CoordinateMetric>) {
    let mut rng = Rng::new(17);
    let data = (0..dimensions)
        .map(|i| match i % 2 {
            0 => rng.normal_vec(FRAMES, 1.0, 0.3),
            _ => rng.uniform_vec(FRAMES, 0.0, TAU),
        })
        .collect();
    let metrics = (0..dimensions)
        .map(|i| match i % 2 {
            0 => CoordinateMetric::Linear,
            _ => CoordinateMetric::Periodic { period: TAU },
        })
        .collect();
    (data, metrics)
}

#[test]
fn mie_results_do_not_depend_on_threads_or_runs() {
    // Enough coordinates that each order's terms (130 singles, 276 pairs, 364
    // triples, 495 quadruples) are summed in several parallel chunks.
    for (order, dimensions) in [(1, 130), (2, 24), (3, 14), (4, 12)] {
        let (data, metrics) = system(dimensions);
        let options = EntropyOptions {
            metrics: Some(&metrics),
            ..Default::default()
        };
        let estimate = || {
            (
                calculate_entropy(&data, FRAMES, order, &options).unwrap(),
                estimate_coordinate_entropy(&data, FRAMES, order, &options).unwrap(),
            )
        };
        let (total, shares) = with_threads(1, estimate);
        for threads in [2, 3, 8, 8, 8] {
            let (other_total, other_shares) = with_threads(threads, estimate);
            assert_eq!(
                total.to_bits(),
                other_total.to_bits(),
                "order {order}, {threads} threads: total {total} vs {other_total}"
            );
            for (a, b) in shares.iter().zip(&other_shares) {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "order {order}, {threads} threads: {a} vs {b}"
                );
            }
        }
    }
}

#[test]
fn network_results_do_not_depend_on_threads_or_runs() {
    let (data, metrics) = system(24);
    let options = EntropyOptions {
        metrics: Some(&metrics),
        ..Default::default()
    };
    let estimate = || {
        (
            estimate_network_entropy(&data, FRAMES, 3, &options).unwrap(),
            estimate_coordinate_network_entropy(&data, FRAMES, 3, &options).unwrap(),
        )
    };
    let (total, shares) = with_threads(1, estimate);
    for threads in [2, 3, 8, 8, 8] {
        let (other_total, other_shares) = with_threads(threads, estimate);
        assert_eq!(total.to_bits(), other_total.to_bits(), "{threads} threads");
        for (a, b) in shares.iter().zip(&other_shares) {
            assert_eq!(a.to_bits(), b.to_bits(), "{threads} threads: {a} vs {b}");
        }
    }
}
