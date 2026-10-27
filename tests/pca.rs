//! Checks of the experimental PCA entropy (`molentropy::pca`).

mod common;

use common::Rng;
use molentropy::pca::{PcaOptions, estimate_pca_entropy};

/// A rigid six-atom reference with independent Gaussian displacements of
/// `sigma` A on every atom, `n` frames.
fn harmonic_frames(rng: &mut Rng, n: usize, sigma: f64) -> Vec<Vec<[f64; 3]>> {
    let reference = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.0, 1.4, 0.0],
        [3.5, 1.5, 0.3],
        [4.1, 2.8, -0.4],
        [5.6, 2.9, 0.1],
    ];
    (0..n)
        .map(|_| {
            reference
                .iter()
                .map(|r| {
                    [
                        r[0] + sigma * rng.normal(),
                        r[1] + sigma * rng.normal(),
                        r[2] + sigma * rng.normal(),
                    ]
                })
                .collect()
        })
        .collect()
}

/// `frame` rotated by a random rotation and shifted by a random translation.
fn moved(rng: &mut Rng, frame: &[[f64; 3]]) -> Vec<[f64; 3]> {
    // A random unit quaternion.
    let q: Vec<f64> = (0..4).map(|_| rng.normal()).collect();
    let norm = q.iter().map(|v| v * v).sum::<f64>().sqrt();
    let (w, x, y, z) = (q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm);
    let r = [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ];
    let shift = [
        10.0 * rng.normal(),
        10.0 * rng.normal(),
        10.0 * rng.normal(),
    ];
    frame
        .iter()
        .map(|p| {
            let mut out = [0.0; 3];
            for i in 0..3 {
                out[i] = r[i][0] * p[0] + r[i][1] * p[1] + r[i][2] * p[2] + shift[i];
            }
            out
        })
        .collect()
}

const MASSES: [f64; 6] = [12.0, 14.0, 12.0, 16.0, 12.0, 1.0];

#[test]
fn rigid_body_motion_does_not_change_the_entropy() {
    let mut rng = Rng::new(41);
    let frames = harmonic_frames(&mut rng, 3_000, 0.3);
    let tumbling: Vec<_> = frames.iter().map(|f| moved(&mut rng, f)).collect();
    let options = PcaOptions::default();
    let still = estimate_pca_entropy(&frames, &MASSES, &options).unwrap();
    let moving = estimate_pca_entropy(&tumbling, &MASSES, &options).unwrap();
    assert!(
        (still.total - moving.total).abs() < 1e-6 * still.total.abs().max(1.0),
        "{} vs {}",
        still.total,
        moving.total
    );
    assert_eq!(still.eigenvalues.len(), 3 * MASSES.len() - 6);
}

#[test]
fn quasi_harmonic_term_is_the_sum_over_the_eigenvalues() {
    let mut rng = Rng::new(42);
    let frames = harmonic_frames(&mut rng, 2_000, 0.2);
    let options = PcaOptions::default();
    let result = estimate_pca_entropy(&frames, &MASSES, &options).unwrap();
    // alpha = hbar / sqrt(k_B T lambda), lambda in amu A^2.
    let thermal = 1.380_649e-23 * options.temperature;
    let expected: f64 = result
        .eigenvalues
        .iter()
        .map(|&lambda| {
            let a = 1.054_571_817e-34 / (thermal * lambda * 1.660_539_066_60e-47).sqrt();
            a / (a.exp() - 1.0) - (1.0 - (-a).exp()).ln()
        })
        .sum();
    assert!(
        (result.quasi_harmonic - expected).abs() < 1e-9,
        "{} vs {expected}",
        result.quasi_harmonic
    );
    assert!(result.eigenvalues.windows(2).all(|w| w[0] >= w[1]));
}

#[test]
fn gaussian_fluctuations_need_little_correction() {
    // Independent Gaussian displacements: after superposition the modes are
    // close to Gaussian and nearly independent, so both corrections are small
    // next to the quasi-harmonic entropy.
    let mut rng = Rng::new(43);
    let frames = harmonic_frames(&mut rng, 5_000, 0.3);
    let result = estimate_pca_entropy(&frames, &MASSES, &PcaOptions::default()).unwrap();
    assert!(result.classical_modes > 0);
    assert!(result.anharmonic >= 0.0 && result.mutual_information >= 0.0);
    let corrections = result.anharmonic + result.mutual_information;
    assert!(
        corrections < 0.05 * result.quasi_harmonic,
        "corrections {corrections} against quasi-harmonic {}",
        result.quasi_harmonic
    );
}

#[test]
fn invalid_input_is_rejected() {
    let mut rng = Rng::new(44);
    let frames = harmonic_frames(&mut rng, 100, 0.3);
    let options = PcaOptions::default();
    assert!(estimate_pca_entropy(&frames, &MASSES[..2], &options).is_err());
    assert!(estimate_pca_entropy(&frames[..5], &MASSES, &options).is_err());
    let mut bad = MASSES;
    bad[0] = 0.0;
    assert!(estimate_pca_entropy(&frames, &bad, &options).is_err());
    let cold = PcaOptions {
        temperature: 0.0,
        ..options
    };
    assert!(estimate_pca_entropy(&frames, &MASSES, &cold).is_err());
}
