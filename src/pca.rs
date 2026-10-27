//! **Experimental.** Entropy from the principal components of the mass-weighted
//! Cartesian coordinates, with corrections for anharmonicity and pairwise
//! dependence (Numata, Wan and Knapp, *Genome Informatics* 18, 192, 2007).
//!
//! This estimator is investigational: it is in the crate so it can be compared
//! with the others on the same benchmarks, not because it is recommended. It
//! works in a different coordinate system from every other estimator here and
//! has not been validated against reference entropies; see the user guide.
//!
//! The method:
//!
//! 1. Remove overall translation and rotation by superposing every frame on
//!    the average structure (mass-weighted, iterated until the average stops
//!    changing).
//! 2. Diagonalize the covariance of the mass-weighted coordinates
//!    `sqrt(m) x`. Of the `3N` eigenvalues `lambda_i` (variances, amu A^2),
//!    the six smallest belong to the removed rigid-body motion and are
//!    dropped.
//! 3. Treat each remaining principal component as a quantum harmonic
//!    oscillator of frequency `omega_i = sqrt(k_B T / lambda_i)`. With
//!    `alpha_i = hbar omega_i / k_B T`, its entropy is
//!    `alpha_i / (e^alpha_i - 1) - ln(1 - e^-alpha_i)`. The sum is the
//!    quasi-harmonic entropy, an upper bound.
//! 4. For the classical modes (`alpha_i < 1`), project the trajectory on them
//!    and subtract
//!    - each mode's anharmonicity: its Gaussian entropy
//!      `1/2 ln(2 pi e lambda_i)` minus its k-nearest-neighbor entropy; a mode
//!      whose correction is below `harmonic_threshold` of its quantum entropy
//!      counts as harmonic (correction 0);
//!    - the mutual information of each pair of modes from one- and
//!      two-dimensional k-nearest-neighbor entropies; negative values count as
//!      0.
//!
//! The k-nearest-neighbor entropy of `n` points in `d` dimensions is
//! `psi(n) - psi(k) + ln V_d + (d / n) sum ln R_k`, with `R_k` the Euclidean
//! distance to the k-th nearest neighbor and `V_d` the volume of the unit
//! ball (Kozachenko and Leonenko; Hnizdo et al.).
//!
//! All results are in units of k_B. Unlike the other estimators, the
//! quasi-harmonic term is an absolute quantum entropy of the vibrations (no
//! momentum part), so only differences between ensembles of the same molecule
//! compare with the rest of the crate.

use std::num::NonZero;

use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::{DMatrix, Matrix3, Vector3};
use rayon::prelude::*;

use crate::log_unit_ball_volume;

/// Reduced Planck constant, J s.
const HBAR: f64 = 1.054_571_817e-34;
/// Boltzmann constant, J/K.
const BOLTZMANN: f64 = 1.380_649e-23;
/// One amu A^2 in kg m^2.
const AMU_ANGSTROM2: f64 = 1.660_539_066_60e-27 * 1e-20;
/// Superposition stops when the average structure moves less than this (A, RMS).
const SUPERPOSITION_TOLERANCE: f64 = 1e-6;
const SUPERPOSITION_ITERATIONS: usize = 50;

/// Settings for [`estimate_pca_entropy`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcaOptions {
    /// Temperature of the simulation, K.
    pub temperature: f64,
    /// The k of the k-nearest-neighbor entropies.
    pub neighbors: usize,
    /// A classical mode whose anharmonic correction is below this fraction of
    /// its quantum entropy counts as harmonic.
    pub harmonic_threshold: f64,
}

impl Default for PcaOptions {
    /// 300 K, the 4th nearest neighbor and a threshold of 0.007, as in the
    /// original method.
    fn default() -> Self {
        Self {
            temperature: 300.0,
            neighbors: 4,
            harmonic_threshold: 0.007,
        }
    }
}

/// The PCA entropy and its parts, in units of k_B.
#[derive(Debug, Clone, PartialEq)]
pub struct PcaEntropy {
    /// `quasi_harmonic - anharmonic - mutual_information`.
    pub total: f64,
    /// Sum of the quantum harmonic-oscillator entropies of the principal
    /// components.
    pub quasi_harmonic: f64,
    /// Sum of the anharmonic corrections, never negative.
    pub anharmonic: f64,
    /// Sum of the pairwise mutual information between classical modes, never
    /// negative.
    pub mutual_information: f64,
    /// Number of classical modes (`alpha < 1`), which the corrections cover.
    pub classical_modes: usize,
    /// Variances of the principal components, amu A^2, largest first; the six
    /// rigid-body ones are left out.
    pub eigenvalues: Vec<f64>,
}

/// The PCA entropy of `frames` (`frames[frame][atom]`, A) of atoms with
/// `masses` (amu). Needs at least three atoms and more frames than the
/// neighbor count.
pub fn estimate_pca_entropy(
    frames: &[Vec<[f64; 3]>],
    masses: &[f64],
    options: &PcaOptions,
) -> Result<PcaEntropy, String> {
    let n_atoms = masses.len();
    if n_atoms < 3 {
        return Err(format!(
            "the PCA entropy needs at least 3 atoms, got {n_atoms}"
        ));
    }
    if masses.iter().any(|&m| !(m > 0.0 && m.is_finite())) {
        return Err("every mass must be positive".to_string());
    }
    if frames.iter().any(|frame| frame.len() != n_atoms) {
        return Err(format!("every frame must have {n_atoms} atoms"));
    }
    if frames.len() <= options.neighbors + 1 {
        return Err(format!(
            "the PCA entropy needs more than {} frames, got {}",
            options.neighbors + 1,
            frames.len()
        ));
    }
    if options.neighbors == 0 {
        return Err("the neighbor count must be at least 1".to_string());
    }
    if !(options.temperature > 0.0 && options.temperature.is_finite()) {
        return Err(format!("invalid temperature {}", options.temperature));
    }

    let aligned = superpose(frames, masses);
    let n_frames = aligned.len();
    let dimension = 3 * n_atoms;
    // Mass-weighted coordinates, centered: frames x 3N.
    let root_masses: Vec<f64> = masses.iter().map(|m| m.sqrt()).collect();
    let mut y = DMatrix::<f64>::from_fn(n_frames, dimension, |f, c| {
        root_masses[c / 3] * aligned[f][c / 3][c % 3]
    });
    for c in 0..dimension {
        let mean = y.column(c).mean();
        y.column_mut(c).add_scalar_mut(-mean);
    }
    let covariance = y.transpose() * &y / (n_frames as f64 - 1.0);
    let eigen = covariance.symmetric_eigen();
    let mut modes: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(i, value)| (value, i))
        .collect();
    modes.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    modes.truncate(dimension.saturating_sub(6));

    let thermal = BOLTZMANN * options.temperature;
    let alpha = |variance: f64| HBAR / (thermal * variance * AMU_ANGSTROM2).sqrt();
    let quantum = |variance: f64| {
        if variance <= 0.0 {
            return 0.0;
        }
        let a = alpha(variance);
        let e = (-a).exp();
        a * e / (1.0 - e) - (-e).ln_1p()
    };
    let quasi_harmonic: f64 = modes.iter().map(|&(variance, _)| quantum(variance)).sum();

    // Projections of the trajectory on the classical modes.
    let classical: Vec<(f64, usize)> = modes
        .iter()
        .copied()
        .filter(|&(variance, _)| variance > 0.0 && alpha(variance) < 1.0)
        .collect();
    let projections: Vec<Vec<f64>> = classical
        .iter()
        .map(|&(_, index)| {
            (&y * eigen.eigenvectors.column(index))
                .iter()
                .copied()
                .collect()
        })
        .collect();
    let k = options.neighbors;
    let marginals: Vec<f64> = projections
        .par_iter()
        .map(|z| knn_entropy::<1>(&z.iter().map(|&v| [v]).collect::<Vec<_>>(), k))
        .collect::<Result<_, _>>()?;
    let anharmonic: f64 = classical
        .iter()
        .zip(&marginals)
        .map(|(&(variance, _), &marginal)| {
            let gaussian = 0.5 * (2.0 * std::f64::consts::PI * std::f64::consts::E * variance).ln();
            let correction = gaussian - marginal;
            if correction < options.harmonic_threshold * quantum(variance) {
                0.0
            } else {
                correction
            }
        })
        .sum();
    let count = classical.len();
    let pairs: Vec<(usize, usize)> = (0..count)
        .flat_map(|i| ((i + 1)..count).map(move |j| (i, j)))
        .collect();
    let pair_information: Vec<f64> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let points: Vec<[f64; 2]> = projections[i]
                .iter()
                .zip(&projections[j])
                .map(|(&a, &b)| [a, b])
                .collect();
            knn_entropy::<2>(&points, k).map(|joint| (marginals[i] + marginals[j] - joint).max(0.0))
        })
        .collect::<Result<_, _>>()?;
    let mutual_information: f64 = pair_information.iter().sum();

    Ok(PcaEntropy {
        total: quasi_harmonic - anharmonic - mutual_information,
        quasi_harmonic,
        anharmonic,
        mutual_information,
        classical_modes: count,
        eigenvalues: modes.iter().map(|&(variance, _)| variance).collect(),
    })
}

/// Every frame superposed on the mass-weighted average structure, centered on
/// its center of mass: fit to the first frame, then refit to the average
/// until it stops moving.
fn superpose(frames: &[Vec<[f64; 3]>], masses: &[f64]) -> Vec<Vec<[f64; 3]>> {
    let total_mass: f64 = masses.iter().sum();
    let centered: Vec<Vec<[f64; 3]>> = frames
        .iter()
        .map(|frame| {
            let mut center = [0.0; 3];
            for (x, &m) in frame.iter().zip(masses) {
                for axis in 0..3 {
                    center[axis] += m * x[axis] / total_mass;
                }
            }
            frame
                .iter()
                .map(|x| [x[0] - center[0], x[1] - center[1], x[2] - center[2]])
                .collect()
        })
        .collect();
    let mut reference = centered[0].clone();
    let mut aligned = centered.clone();
    for _ in 0..SUPERPOSITION_ITERATIONS {
        aligned = centered
            .par_iter()
            .map(|frame| rotate_onto(frame, &reference, masses))
            .collect();
        let n = aligned.len() as f64;
        let average: Vec<[f64; 3]> = (0..masses.len())
            .map(|atom| {
                let mut mean = [0.0; 3];
                for frame in &aligned {
                    for axis in 0..3 {
                        mean[axis] += frame[atom][axis] / n;
                    }
                }
                mean
            })
            .collect();
        let shift = (average
            .iter()
            .zip(&reference)
            .map(|(a, b)| (0..3).map(|axis| (a[axis] - b[axis]).powi(2)).sum::<f64>())
            .sum::<f64>()
            / masses.len() as f64)
            .sqrt();
        reference = average;
        if shift < SUPERPOSITION_TOLERANCE {
            break;
        }
    }
    aligned
}

/// `frame` (centered) rotated to minimize its mass-weighted squared distance
/// to `reference` (centered), by the Kabsch method.
fn rotate_onto(frame: &[[f64; 3]], reference: &[[f64; 3]], masses: &[f64]) -> Vec<[f64; 3]> {
    let mut h = Matrix3::<f64>::zeros();
    for ((x, r), &m) in frame.iter().zip(reference).zip(masses) {
        h += m * Vector3::from(*x) * Vector3::from(*r).transpose();
    }
    let svd = h.svd(true, true);
    let (u, v_t) = (svd.u.unwrap(), svd.v_t.unwrap());
    let mut correction = Matrix3::<f64>::identity();
    if (v_t.transpose() * u.transpose()).determinant() < 0.0 {
        correction[(2, 2)] = -1.0;
    }
    let rotation = v_t.transpose() * correction * u.transpose();
    frame
        .iter()
        .map(|x| {
            let rotated = rotation * Vector3::from(*x);
            [rotated[0], rotated[1], rotated[2]]
        })
        .collect()
}

/// k-nearest-neighbor entropy of `points` in `D` dimensions, Euclidean.
fn knn_entropy<const D: usize>(points: &[[f64; D]], k: usize) -> Result<f64, String> {
    let n = points.len();
    let tree: ImmutableKdTree<f64, D> = ImmutableKdTree::new_from_slice(points);
    let log_distances: Vec<f64> = points
        .par_iter()
        .map(|point| {
            // The point itself comes first, at distance 0.
            let neighbors =
                tree.nearest_n::<SquaredEuclidean>(point, NonZero::new(k + 1).unwrap());
            let squared = neighbors.last().map_or(0.0, |neighbor| neighbor.distance);
            if squared > 0.0 {
                Ok(0.5 * squared.ln())
            } else {
                Err(format!(
                    "a point has {k} or more identical copies; the {k}-nearest-neighbor entropy is undefined"
                ))
            }
        })
        .collect::<Result<_, _>>()?;
    let mean: f64 = log_distances.iter().sum::<f64>() / n as f64;
    // psi(n) - psi(k) = H_(n-1) - H_(k-1).
    let digamma_difference: f64 = (k..n).map(|i| 1.0 / i as f64).sum();
    Ok(digamma_difference + log_unit_ball_volume(D) + D as f64 * mean)
}

#[cfg(test)]
mod tests {
    use super::knn_entropy;

    /// Normal samples by the Box-Muller transform from a fixed linear
    /// congruential sequence.
    fn normals(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        let mut uniform = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        };
        (0..n)
            .map(|_| (-2.0 * uniform().ln()).sqrt() * (std::f64::consts::TAU * uniform()).cos())
            .collect()
    }

    #[test]
    fn knn_entropy_of_a_gaussian_matches_the_closed_form() {
        let exact = |d: f64| 0.5 * d * (2.0 * std::f64::consts::PI * std::f64::consts::E).ln();
        let x = normals(20_000, 1);
        let y = normals(20_000, 2);
        for k in [1, 4] {
            let one: Vec<[f64; 1]> = x.iter().map(|&v| [v]).collect();
            let two: Vec<[f64; 2]> = x.iter().zip(&y).map(|(&a, &b)| [a, b]).collect();
            let e1 = knn_entropy(&one, k).unwrap();
            let e2 = knn_entropy(&two, k).unwrap();
            assert!(
                (e1 - exact(1.0)).abs() < 0.02,
                "k {k}: 1D {e1} vs {}",
                exact(1.0)
            );
            assert!(
                (e2 - exact(2.0)).abs() < 0.03,
                "k {k}: 2D {e2} vs {}",
                exact(2.0)
            );
        }
    }
}
