//! Theiler-window exclusion of temporally adjacent neighbors.
//!
//! A window `w` restricts frame i's nearest neighbor to frames j with
//! |i - j| > w, and replaces the sample term psi(N) with the mean over frames
//! of psi(n_i), where n_i is N minus the other frames in frame i's window.
//! These tests check that `w = 0` is the unwindowed estimator bit for bit, that
//! the kd-tree search agrees with a brute-force reference for every window,
//! that the constant has the closed form above, and that the window removes the
//! spurious correlation an oversampled series feeds the estimator.

mod common;

use common::{Rng, brute_force_log_nn_sum_windowed, normal_entropy};
use molentropy::{
    CoordinateMetric, EntropyOptions, FiniteSampleConstant, calc_joint_nn,
    calc_joint_nn_with_window, calculate_entropy, estimate_coordinate_entropy,
    estimate_coordinate_mutual_information,
};

const EULER_MASCHERONI: f64 = std::f64::consts::EULER_GAMMA;
const TAU: f64 = std::f64::consts::TAU;

fn windowed(theiler_window: usize, constant: FiniteSampleConstant) -> EntropyOptions<'static> {
    EntropyOptions {
        constant,
        theiler_window,
        ..Default::default()
    }
}

/// Random series quantized to a coarse grid so that exact duplicates, which
/// the estimator must skip, occur both inside and outside the window.
fn quantized_series(rng: &mut Rng, n: usize, high: f64, step: f64) -> Vec<f64> {
    (0..n)
        .map(|_| (rng.uniform_range(0.0, high) / step).round() * step)
        .collect()
}

fn assert_close(label: &str, actual: f64, expected: f64) {
    let tolerance = 1e-9 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() < tolerance,
        "{label}: {actual} vs {expected}"
    );
}

// ---------------------------------------------------------------------------
// w = 0 is the unwindowed estimator
// ---------------------------------------------------------------------------

#[test]
fn zero_window_is_bitwise_identical_to_unwindowed_search() {
    let mut rng = Rng::new(11);
    let a = quantized_series(&mut rng, 300, TAU, 0.05);
    let b = rng.normal_vec(300, 0.0, 1.0);
    let c = quantized_series(&mut rng, 300, 3.0, 0.1);
    let d = rng.uniform_vec(300, 0.0, TAU);
    let linear = CoordinateMetric::Linear;
    let periodic = CoordinateMetric::Periodic { period: TAU };

    assert_eq!(
        calc_joint_nn([&a], [periodic]).unwrap().to_bits(),
        calc_joint_nn_with_window([&a], [periodic], 0)
            .unwrap()
            .to_bits()
    );
    assert_eq!(
        calc_joint_nn([&a, &b], [periodic, linear])
            .unwrap()
            .to_bits(),
        calc_joint_nn_with_window([&a, &b], [periodic, linear], 0)
            .unwrap()
            .to_bits()
    );
    assert_eq!(
        calc_joint_nn([&a, &b, &c, &d], [periodic, linear, linear, periodic])
            .unwrap()
            .to_bits(),
        calc_joint_nn_with_window([&a, &b, &c, &d], [periodic, linear, linear, periodic], 0)
            .unwrap()
            .to_bits()
    );
}

#[test]
fn explicit_zero_window_matches_default_options() {
    let mut rng = Rng::new(12);
    let data: Vec<Vec<f64>> = (0..5).map(|_| rng.normal_vec(400, 0.0, 1.0)).collect();
    for constant in [
        FiniteSampleConstant::Asymptotic,
        FiniteSampleConstant::Exact,
    ] {
        let default = EntropyOptions {
            constant,
            ..Default::default()
        };
        for order in 1..=4 {
            assert_eq!(
                calculate_entropy(&data, 400, order, &default)
                    .unwrap()
                    .to_bits(),
                calculate_entropy(&data, 400, order, &windowed(0, constant))
                    .unwrap()
                    .to_bits(),
                "{constant:?} order {order}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The kd-tree search against brute force
// ---------------------------------------------------------------------------

#[test]
fn windowed_search_matches_brute_force_in_every_dimension() {
    let linear = CoordinateMetric::Linear;
    let periodic = CoordinateMetric::Periodic { period: TAU };
    for (seed, window) in [(21, 1), (22, 2), (23, 7), (24, 40)] {
        let mut rng = Rng::new(seed);
        // An AR(1) series makes temporal neighbors also spatial neighbors, the
        // case where the window actually changes which point is nearest.
        let a: Vec<f64> = rng
            .ar1_vec(250, 0.95)
            .into_iter()
            .map(|x| x * 2.0)
            .collect();
        let b = quantized_series(&mut rng, 250, TAU, 0.2);
        let c = rng.ar1_vec(250, 0.8);
        let d = rng.uniform_vec(250, -1.0, 1.0);

        let check = |label: &str, actual: f64, coords: &[&[f64]], metrics: &[CoordinateMetric]| {
            let expected = brute_force_log_nn_sum_windowed(coords, metrics, window);
            assert_close(&format!("{label} w={window}"), actual, expected);
        };
        check(
            "1D linear",
            calc_joint_nn_with_window([&a], [linear], window).unwrap(),
            &[&a],
            &[linear],
        );
        check(
            "1D periodic",
            calc_joint_nn_with_window([&b], [periodic], window).unwrap(),
            &[&b],
            &[periodic],
        );
        check(
            "2D mixed",
            calc_joint_nn_with_window([&a, &b], [linear, periodic], window).unwrap(),
            &[&a, &b],
            &[linear, periodic],
        );
        check(
            "3D mixed",
            calc_joint_nn_with_window([&b, &c, &d], [periodic, linear, linear], window).unwrap(),
            &[&b, &c, &d],
            &[periodic, linear, linear],
        );
        check(
            "4D mixed",
            calc_joint_nn_with_window([&a, &b, &c, &d], [linear, periodic, linear, linear], window)
                .unwrap(),
            &[&a, &b, &c, &d],
            &[linear, periodic, linear, linear],
        );
    }
}

#[test]
fn window_skips_duplicates_at_every_distance() {
    // Frames 0 and 3 coincide. With w = 1 frame 0 may only use frames 2..5:
    // frame 3 is a duplicate, so its nearest distinct eligible point is frame 2
    // at distance 1, not frame 1 at distance 0.5 (inside the window).
    let series = [0.0, 0.5, 1.0, 0.0, 4.0, 6.0];
    let actual = calc_joint_nn_with_window([&series], [CoordinateMetric::Linear], 1).unwrap();
    let expected = brute_force_log_nn_sum_windowed(&[&series], &[CoordinateMetric::Linear], 1);
    assert_close("duplicates", actual, expected);
    // Frame 0 -> frame 2 (1.0); frame 1 -> frame 3 (0.5); frame 2 -> frame 0
    // (1.0); frame 3 -> frame 1 (0.5); frame 4 -> frame 2 (3.0); frame 5 ->
    // frame 2 (5.0).
    let by_hand = [1.0_f64, 0.5, 1.0, 0.5, 3.0, 5.0]
        .iter()
        .map(|d| d.ln())
        .sum();
    assert_close("duplicates by hand", actual, by_hand);
}

// ---------------------------------------------------------------------------
// Closed-form constant
// ---------------------------------------------------------------------------

#[test]
fn windowed_constant_matches_closed_form() {
    // Four frames one unit apart with w = 1. Every frame's nearest eligible
    // neighbor is two units away, so the mean log distance is ln 2. Frames 0
    // and 3 lose one neighbor to the window (n_i = 3), frames 1 and 2 lose two
    // (n_i = 2).
    let data = vec![vec![0.0, 1.0, 2.0, 3.0]];
    let log_volume = 2.0_f64.ln();

    let exact = calculate_entropy(&data, 4, 1, &windowed(1, FiniteSampleConstant::Exact)).unwrap();
    // Mean of H_(n_i - 1): (H_2 + H_1 + H_1 + H_2) / 4 = (1.5 + 1 + 1 + 1.5) / 4.
    assert_close("exact", exact, 2.0_f64.ln() + log_volume + 1.25);

    let asymptotic =
        calculate_entropy(&data, 4, 1, &windowed(1, FiniteSampleConstant::Asymptotic)).unwrap();
    // Mean of ln(n_i) + gamma: (2 ln 3 + 2 ln 2) / 4 + gamma.
    assert_close(
        "asymptotic",
        asymptotic,
        2.0_f64.ln() + log_volume + 6.0_f64.ln() / 2.0 + EULER_MASCHERONI,
    );
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

#[test]
fn window_needs_two_w_plus_two_frames() {
    let ok = vec![vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]];
    calculate_entropy(&ok, 6, 1, &windowed(2, FiniteSampleConstant::Exact))
        .expect("six frames are enough for w = 2");

    let error = calculate_entropy(&ok, 5, 1, &windowed(2, FiniteSampleConstant::Exact))
        .expect_err("five frames leave frame 2 without an eligible neighbor");
    assert!(
        error.contains("Theiler window 2 needs at least 6 frames"),
        "{error}"
    );

    let error = calc_joint_nn_with_window([&ok[0]], [CoordinateMetric::Linear], 3)
        .expect_err("w = 3 needs eight frames");
    assert!(error.contains("Theiler window 3"), "{error}");
}

#[test]
fn window_errors_when_only_duplicates_lie_outside_it() {
    // Frame 0's only eligible frames (2 and 3) sit at its own position.
    let series = [0.0, 5.0, 0.0, 0.0];
    let error = calc_joint_nn_with_window([&series], [CoordinateMetric::Linear], 1)
        .expect_err("frame 0 has no distinct eligible neighbor");
    assert!(error.contains("frame 0"), "{error}");
}

#[test]
fn every_estimator_honors_the_window() {
    // The per-coordinate estimators must sum to the total under the window, as
    // they do without it.
    let mut rng = Rng::new(31);
    let data: Vec<Vec<f64>> = (0..4).map(|_| rng.ar1_vec(2000, 0.95)).collect();
    let options = windowed(10, FiniteSampleConstant::Exact);
    for order in 1..=2 {
        let total = calculate_entropy(&data, 2000, order, &options).unwrap();
        let shares: f64 = estimate_coordinate_entropy(&data, 2000, order, &options)
            .unwrap()
            .iter()
            .sum();
        assert_close(&format!("order {order} shares"), shares, total);
    }
    let unwindowed = calculate_entropy(&data, 2000, 1, &windowed(0, FiniteSampleConstant::Exact));
    assert_ne!(
        unwindowed.unwrap(),
        calculate_entropy(&data, 2000, 1, &options).unwrap()
    );
}

// ---------------------------------------------------------------------------
// Statistical behavior
// ---------------------------------------------------------------------------

/// Mean over replicates of `estimate(seed)`.
fn replicate_mean(replicates: u64, estimate: impl Fn(u64) -> f64) -> f64 {
    (0..replicates).map(&estimate).sum::<f64>() / replicates as f64
}

#[test]
fn window_leaves_independent_samples_unbiased() {
    // With iid samples the window only discards candidates, and the per-frame
    // constant accounts for that, so the estimate must not move systematically.
    let options = |w| windowed(w, FiniteSampleConstant::Exact);
    for window in [50, 500] {
        let shift = replicate_mean(8, |seed| {
            let mut rng = Rng::new(4000 + seed);
            let data: Vec<Vec<f64>> = (0..2).map(|_| rng.normal_vec(10_000, 0.0, 1.0)).collect();
            calculate_entropy(&data, 10_000, 2, &options(window)).unwrap()
                - calculate_entropy(&data, 10_000, 2, &options(0)).unwrap()
        });
        println!("iid order-2 shift at w={window}: {shift:+.5}");
        assert!(
            shift.abs() < 0.01,
            "w={window} shifted iid entropy by {shift}"
        );
    }
}

#[test]
fn window_removes_spurious_mutual_information_of_oversampled_series() {
    // Two independent AR(1) series with lag-one correlation 0.99 have zero
    // mutual information and standard normal marginals. Without a window,
    // each frame's nearest neighbor is its own predecessor or successor in
    // both coordinates at once, which reads as dependence. Measured over 32
    // replicates: MI +0.044 +- 0.004 at w = 0 and -0.002 +- 0.004 at w = 5.
    //
    // A larger window over-corrects on this strongly oversampled series (MI
    // about -0.04 at w >= 50, from the clustering of the remaining
    // candidates), so this asserts the bias reduction near the decorrelation
    // scale rather than convergence as w grows.
    let mi = |window| {
        replicate_mean(8, |seed| {
            let mut rng = Rng::new(7000 + seed);
            let data: Vec<Vec<f64>> = (0..2).map(|_| rng.ar1_vec(20_000, 0.99)).collect();
            estimate_coordinate_mutual_information(
                &data,
                20_000,
                &windowed(window, FiniteSampleConstant::Exact),
            )
            .unwrap()[0]
        })
    };
    let unwindowed = mi(0);
    let windowed_mi = mi(5);
    println!("AR(1) phi=0.99 MI: w=0 {unwindowed:+.4}, w=5 {windowed_mi:+.4}");
    assert!(
        unwindowed > 0.03,
        "expected spurious MI without a window, got {unwindowed}"
    );
    assert!(windowed_mi.abs() < 0.015, "window left MI {windowed_mi}");
}

#[test]
fn window_does_not_disturb_marginal_entropy_of_mildly_correlated_series() {
    // phi = 0.9: the marginal estimate is already close at w = 0 and must stay
    // close to 1/2 ln(2 pi e) once the window is applied.
    let estimate = replicate_mean(8, |seed| {
        let mut rng = Rng::new(5000 + seed);
        let data = vec![rng.ar1_vec(20_000, 0.9)];
        calculate_entropy(&data, 20_000, 1, &windowed(20, FiniteSampleConstant::Exact)).unwrap()
    });
    let error = estimate - normal_entropy(1.0);
    println!("AR(1) phi=0.9 H at w=20: error {error:+.4}");
    assert!(error.abs() < 0.01, "error {error}");
}
