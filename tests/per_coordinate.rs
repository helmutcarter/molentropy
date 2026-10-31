//! Per-coordinate shares of the MIE at orders 3 and 4, of MIST and of the
//! split estimator: each set of shares sums to its total, follows the
//! coordinates when they are reordered, and matches an independent
//! computation or a closed form.

mod common;

use std::f64::consts::TAU;

use common::Rng;
use molentropy::quasi_harmonic::{PeriodicCorrelation, QuasiHarmonicOptions};
use molentropy::split::{
    Coupling, HardEstimator, SoftEstimator, SoftFeatures, SoftSelection, SplitOptions,
    StateDefinition, StateOptions, estimate_coordinate_split_entropy, estimate_split_entropy,
};
use molentropy::{
    CoordinateMetric, EntropyOptions, calculate_entropy, estimate_coordinate_entropy,
    estimate_coordinate_mist_entropy, estimate_coordinate_mutual_information,
    estimate_coordinate_network_entropy, estimate_joint_entropy, estimate_mist_entropy,
    estimate_network_entropy,
};

const TORSION: CoordinateMetric = CoordinateMetric::Periodic { period: TAU };
const LINEAR: CoordinateMetric = CoordinateMetric::Linear;

fn options(metrics: &[CoordinateMetric]) -> EntropyOptions<'_> {
    EntropyOptions {
        metrics: Some(metrics),
        ..Default::default()
    }
}

/// Five coordinates: two correlated normals, a torsion that one of them
/// follows, an independent torsion, and an independent normal.
fn five(rng: &mut Rng, n: usize) -> (Vec<Vec<f64>>, Vec<CoordinateMetric>) {
    let (a, b) = rng.correlated_normal_pair(n, 0.7);
    let phi: Vec<f64> = a
        .iter()
        .map(|x| (1.5 * x + 0.5 * rng.normal()).rem_euclid(TAU))
        .collect();
    let free = rng.uniform_vec(n, 0.0, TAU);
    let alone = rng.normal_vec(n, 0.0, 2.0);
    (
        vec![a, b, phi, free, alone],
        vec![LINEAR, LINEAR, TORSION, TORSION, LINEAR],
    )
}

fn close(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs()))
}

// ---------------------------------------------------------------------------
// MIE at orders 3 and 4
// ---------------------------------------------------------------------------

#[test]
fn higher_order_shares_sum_to_the_total_and_follow_the_coordinates() {
    let mut rng = Rng::new(301);
    let (data, metrics) = five(&mut rng, 1500);
    for order in [3, 4] {
        let shares = estimate_coordinate_entropy(&data, 1500, order, &options(&metrics)).unwrap();
        let total = calculate_entropy(&data, 1500, order, &options(&metrics)).unwrap();
        assert!(
            close(shares.iter().sum(), total, 1e-10),
            "order {order}: {} vs {total}",
            shares.iter().sum::<f64>()
        );

        let permutation = [3, 0, 4, 2, 1];
        let permuted: Vec<Vec<f64>> = permutation.iter().map(|&i| data[i].clone()).collect();
        let permuted_metrics: Vec<CoordinateMetric> =
            permutation.iter().map(|&i| metrics[i]).collect();
        let permuted_shares =
            estimate_coordinate_entropy(&permuted, 1500, order, &options(&permuted_metrics))
                .unwrap();
        for (position, &i) in permutation.iter().enumerate() {
            assert!(
                close(permuted_shares[position], shares[i], 1e-10),
                "order {order}, coordinate {i}"
            );
        }
    }
}

#[test]
fn order_three_splits_each_interaction_term_evenly() {
    // An independent computation: the interaction terms m(S) by inclusion-
    // exclusion over joint entropies of subsets, each split evenly among the
    // coordinates of S. With three coordinates at order 3 every subset counts.
    let mut rng = Rng::new(302);
    let (data, metrics) = five(&mut rng, 1200);
    let (data, metrics) = (data[..3].to_vec(), metrics[..3].to_vec());
    let joint = |subset: &[usize]| {
        let columns: Vec<Vec<f64>> = subset.iter().map(|&i| data[i].clone()).collect();
        let subset_metrics: Vec<CoordinateMetric> = subset.iter().map(|&i| metrics[i]).collect();
        estimate_joint_entropy(&columns, 1200, &options(&subset_metrics)).unwrap()
    };
    let h: Vec<f64> = (0..3).map(|i| joint(&[i])).collect();
    let pair = |i: usize, j: usize| joint(&[i, j]) - h[i] - h[j];
    let triple =
        joint(&[0, 1, 2]) - joint(&[0, 1]) - joint(&[0, 2]) - joint(&[1, 2]) + h[0] + h[1] + h[2];
    let expected = [
        h[0] + 0.5 * (pair(0, 1) + pair(0, 2)) + triple / 3.0,
        h[1] + 0.5 * (pair(0, 1) + pair(1, 2)) + triple / 3.0,
        h[2] + 0.5 * (pair(0, 2) + pair(1, 2)) + triple / 3.0,
    ];
    let shares = estimate_coordinate_entropy(&data, 1200, 3, &options(&metrics)).unwrap();
    for i in 0..3 {
        assert!(
            close(shares[i], expected[i], 1e-10),
            "{i}: {} vs {}",
            shares[i],
            expected[i]
        );
    }
}

#[test]
fn an_independent_coordinate_keeps_its_own_entropy_at_every_order() {
    // Every interaction term with an independent coordinate vanishes, up to
    // sampling noise, so its share is its marginal entropy. The pairwise
    // estimates are biased slightly upward at small samples (0.065 nats summed
    // over four pairs at 4,000 frames), so this uses 20,000.
    let mut rng = Rng::new(303);
    let (data, metrics) = five(&mut rng, 20_000);
    let marginal = estimate_coordinate_entropy(&data, 20_000, 1, &options(&metrics)).unwrap()[4];
    for order in [2, 3, 4] {
        let share =
            estimate_coordinate_entropy(&data, 20_000, order, &options(&metrics)).unwrap()[4];
        assert!(
            (share - marginal).abs() < 0.05,
            "order {order}: {share} vs {marginal}"
        );
    }
}

#[test]
fn unsupported_orders_are_refused() {
    let data = vec![vec![0.1, 0.4, 0.8], vec![1.0, 1.3, 1.9]];
    for order in [0, 5] {
        let error = estimate_coordinate_entropy(&data, 3, order, &EntropyOptions::default())
            .expect_err("expected rejection");
        assert!(error.contains("unsupported MIE order"), "{error}");
    }
}

// ---------------------------------------------------------------------------
// MIST
// ---------------------------------------------------------------------------

#[test]
fn mist_shares_take_half_of_each_tree_link() {
    // x1 -> x2 -> x3: the tree is the chain, so x2 loses half of both links
    // and the ends half of one each.
    let mut rng = Rng::new(304);
    let n = 3000;
    let x1 = rng.normal_vec(n, 0.0, 1.0);
    let x2: Vec<f64> = x1.iter().map(|x| 0.8 * x + 0.6 * rng.normal()).collect();
    let x3: Vec<f64> = x2.iter().map(|x| 0.8 * x + 0.6 * rng.normal()).collect();
    let data = vec![x1, x2, x3];
    let metrics = [LINEAR; 3];
    let marginal = estimate_coordinate_entropy(&data, n, 1, &options(&metrics)).unwrap();
    let mi = estimate_coordinate_mutual_information(&data, n, &options(&metrics)).unwrap();
    let (i12, i23) = (mi[0], mi[2]);
    assert!(mi[1] < i12.min(i23), "the chain should be the tree");
    let shares = estimate_coordinate_mist_entropy(&data, n, &options(&metrics)).unwrap();
    let expected = [
        marginal[0] - 0.5 * i12,
        marginal[1] - 0.5 * (i12 + i23),
        marginal[2] - 0.5 * i23,
    ];
    for i in 0..3 {
        assert!(close(shares[i], expected[i], 1e-12), "{i}");
    }
    let total = estimate_mist_entropy(&data, n, &options(&metrics)).unwrap();
    assert!(close(shares.iter().sum(), total, 1e-12));
}

#[test]
fn network_shares_sum_to_the_network_entropy() {
    // A collider, x3 = x1 + x2 + noise, so the third coordinate takes two parents.
    let mut rng = Rng::new(305);
    let n = 3000;
    let x1 = rng.normal_vec(n, 0.0, 1.0);
    let x2 = rng.normal_vec(n, 0.0, 1.0);
    let x3: Vec<f64> = x1
        .iter()
        .zip(&x2)
        .map(|(a, b)| 0.6 * (a + b) + 0.5 * rng.normal())
        .collect();
    let data = vec![x1, x2, x3];
    let metrics = [LINEAR; 3];
    let one = estimate_coordinate_network_entropy(&data, n, 1, &options(&metrics)).unwrap();
    let mist = estimate_coordinate_mist_entropy(&data, n, &options(&metrics)).unwrap();
    for i in 0..3 {
        assert!(
            close(one[i], mist[i], 1e-12),
            "{i}: {} vs {}",
            one[i],
            mist[i]
        );
    }
    let shares = estimate_coordinate_network_entropy(&data, n, 2, &options(&metrics)).unwrap();
    let total = estimate_network_entropy(&data, n, 2, &options(&metrics)).unwrap();
    assert!(close(shares.iter().sum(), total, 1e-12));
    assert!(total < estimate_mist_entropy(&data, n, &options(&metrics)).unwrap());
}

// ---------------------------------------------------------------------------
// Split
// ---------------------------------------------------------------------------

/// A torsion with three wells, a hard coordinate whose spread depends on the
/// well, an independent hard torsion, and an independent linear coordinate.
fn split_system(rng: &mut Rng, n: usize) -> (Vec<Vec<f64>>, Vec<CoordinateMetric>, Vec<u32>) {
    let wells = [
        (0.0, 0.5, 0.1),
        (TAU / 3.0, 0.3, 0.2),
        (2.0 * TAU / 3.0, 0.2, 0.4),
    ];
    let (mut phi, mut x, mut labels) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..n {
        let mut pick = rng.uniform();
        let well = wells
            .iter()
            .position(|w| {
                pick -= w.1;
                pick < 0.0
            })
            .unwrap_or(2);
        phi.push((wells[well].0 + 0.2 * rng.normal()).rem_euclid(TAU));
        x.push(0.3 * phi.last().unwrap().cos() + wells[well].2 * rng.normal());
        labels.push(well as u32);
    }
    let narrow: Vec<f64> = (0..n)
        .map(|_| (1.0 + 0.15 * rng.normal()).rem_euclid(TAU))
        .collect();
    let free = rng.normal_vec(n, 0.0, 1.0);
    (
        vec![phi, x, narrow, free],
        vec![TORSION, LINEAR, TORSION, LINEAR],
        labels,
    )
}

#[test]
fn split_shares_sum_to_each_term_for_every_supported_setting() {
    let mut rng = Rng::new(305);
    let (data, metrics, labels) = split_system(&mut rng, 4000);
    let soft = [0];
    let mut settings = Vec::new();
    let base = SplitOptions::new(SoftSelection::Indices(&soft));
    settings.push(base);
    let mut mist = base;
    mist.soft_estimator = SoftEstimator::MutualInformationSpanningTree;
    mist.coupling = Coupling::Gaussian {
        features: SoftFeatures::Complex,
    };
    settings.push(mist);
    let mut bqh = base;
    bqh.hard_estimator = HardEstimator::BoltzmannQuasiHarmonic;
    bqh.hard = QuasiHarmonicOptions {
        periodic: PeriodicCorrelation::Complex,
        ..Default::default()
    };
    settings.push(bqh);
    let mut states = base;
    states.coupling = Coupling::PerState(StateOptions {
        definition: StateDefinition::Labels(&labels),
        ..Default::default()
    });
    settings.push(states);
    let mut detected_states = bqh;
    detected_states.soft = SoftSelection::Indices(&[0, 2]);
    detected_states.soft_estimator = SoftEstimator::NearestNeighborMie { order: 3 };
    detected_states.coupling = Coupling::PerState(StateOptions::default());
    settings.push(detected_states);
    let mut none = base;
    none.coupling = Coupling::None;
    settings.push(none);
    let mut pairwise = base;
    pairwise.hard_estimator = HardEstimator::QuasiHarmonicMie2;
    settings.push(pairwise);
    let mut first_order = base;
    first_order.hard_estimator = HardEstimator::QuasiHarmonicMie1;
    settings.push(first_order);
    let mut pairwise_states = states;
    pairwise_states.hard_estimator = HardEstimator::QuasiHarmonicMie2;
    settings.push(pairwise_states);

    for setting in &settings {
        let shares =
            estimate_coordinate_split_entropy(&data, 4000, setting, &options(&metrics)).unwrap();
        let plain = estimate_split_entropy(&data, 4000, setting, &options(&metrics)).unwrap();
        assert_eq!(shares.split, plain, "{setting:?}");
        let sum = |v: &[f64]| v.iter().sum::<f64>();
        assert!(
            close(sum(&shares.own), plain.soft + plain.hard, 1e-10),
            "{setting:?}"
        );
        assert!(
            close(sum(&shares.coupling), plain.coupling, 1e-10),
            "{setting:?}"
        );
        assert!(close(sum(&shares.total), plain.total, 1e-10), "{setting:?}");
        for i in 0..4 {
            assert_eq!(shares.total[i], shares.own[i] + shares.coupling[i]);
        }
        if setting.coupling == Coupling::None {
            assert!(shares.coupling.iter().all(|&c| c == 0.0));
        }
        if matches!(setting.coupling, Coupling::PerState(_)) {
            for &i in &plain.soft_coordinates {
                assert_eq!(
                    shares.coupling[i], 0.0,
                    "a soft coordinate's per-state share"
                );
            }
        }
    }
}

#[test]
fn gaussian_coupling_of_one_pair_is_shared_evenly() {
    // A free rotor phi and x = a cos(phi) + noise: the coupling is
    // 1/2 ln(1 - rho^2), all of it in the (cos phi, x) block, whose log
    // diagonal gives each of the two 1/4 ln(1 - rho^2). The other columns
    // (sin phi and the higher harmonics) are uncorrelated with x.
    let mut rng = Rng::new(306);
    let n = 50_000;
    let (a, sigma) = (1.0, 0.5);
    let phi = rng.uniform_vec(n, 0.0, TAU);
    let x: Vec<f64> = phi
        .iter()
        .map(|p| a * p.cos() + sigma * rng.normal())
        .collect();
    let data = vec![phi, x];
    let metrics = [TORSION, LINEAR];
    let rho_squared = (a * a / 2.0) / (a * a / 2.0 + sigma * sigma);
    let quarter = 0.25 * (1.0 - rho_squared).ln();
    let shares = estimate_coordinate_split_entropy(
        &data,
        n,
        &SplitOptions::new(SoftSelection::Indices(&[0])),
        &options(&metrics),
    )
    .unwrap();
    for i in 0..2 {
        assert!(
            (shares.coupling[i] - quarter).abs() < 0.01,
            "{i}: {} vs {quarter}",
            shares.coupling[i]
        );
    }
}

#[test]
fn per_state_coupling_falls_on_the_coordinate_that_depends_on_the_state() {
    let mut rng = Rng::new(307);
    let (data, metrics, labels) = split_system(&mut rng, 6000);
    let mut settings = SplitOptions::new(SoftSelection::Indices(&[0]));
    settings.coupling = Coupling::PerState(StateOptions {
        definition: StateDefinition::Labels(&labels),
        ..Default::default()
    });
    let shares =
        estimate_coordinate_split_entropy(&data, 6000, &settings, &options(&metrics)).unwrap();
    // x's spread depends on the well; the narrow torsion and the free normal
    // do not.
    assert!(shares.coupling[1] < -0.2, "{}", shares.coupling[1]);
    for i in [2, 3] {
        assert!(
            shares.coupling[i].abs() < 0.03,
            "{i}: {}",
            shares.coupling[i]
        );
    }
}

#[test]
fn the_joint_soft_estimate_has_no_shares() {
    let mut rng = Rng::new(308);
    let (data, metrics, _) = split_system(&mut rng, 1000);
    let mut settings = SplitOptions::new(SoftSelection::Indices(&[0, 2]));
    settings.soft_estimator = SoftEstimator::JointNearestNeighbor;
    let error =
        estimate_coordinate_split_entropy(&data, 1000, &settings, &options(&metrics)).unwrap_err();
    assert!(error.contains("joint soft estimate"), "{error}");
}
