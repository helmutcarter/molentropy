//! Quasi-harmonic (QH) and Boltzmann quasi-harmonic (BQH) entropy.
//!
//! Both treat correlations through the correlation matrix C of the
//! coordinates, C_ij = sigma_ij / sqrt(sigma_ii sigma_jj), which contributes
//! 1/2 ln|C| <= 0. They differ in the first-order term:
//!
//! - QH (Karplus and Kushick, 1981) fits a Gaussian to each coordinate, which
//!   contributes 1/2 ln(2 pi e sigma_ii).
//! - BQH (Harpole and Sharp, J. Phys. Chem. B 2011) replaces the Gaussian
//!   marginals by the entropy of each coordinate's actual distribution,
//!   -integral p ln p, which matters for multimodal coordinates such as
//!   torsions with several rotamers.
//!
//! All values are differential entropies in units of k_B, like the
//! nearest-neighbor estimators, so they can be compared directly.
//!
//! Periodic coordinates enter C as z = e^(i 2 pi q / P), following Harpole
//! and Sharp, which makes C Hermitian; see [`PeriodicCorrelation`]. For the
//! first-order Gaussian term they are centered on their circular mean and
//! wrapped to (-P/2, P/2], so sigma_ii is in the coordinate's own units.

use rayon::prelude::*;

use crate::{CoordinateMetric, EntropyOptions, coordinate_entropy_impl, resolve_input};

/// How the first-order entropy -integral p ln p of each coordinate is estimated
/// for BQH.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MarginalEstimator {
    /// The one-dimensional nearest-neighbor estimate, the same first-order term
    /// as the MIE. Uses `EntropyOptions::constant` and `theiler_window`.
    #[default]
    NearestNeighbor,
    /// A histogram of `bins` equal bins, -sum p_k ln(p_k / width), Harpole and
    /// Sharp's eq 9. Linear coordinates are binned over their observed range,
    /// periodic ones over a full period. Some later BQH work reports the
    /// discrete entropy -sum p_k ln p_k, which is this value minus ln(width).
    Histogram {
        /// Number of equal-width bins per coordinate.
        bins: usize,
    },
}

/// How periodic coordinates enter the correlation matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PeriodicCorrelation {
    /// As the complex number z = e^(i 2 pi q / P), so C is Hermitian (Harpole
    /// and Sharp). Correlations on the circle, such as the fixed 120 degree
    /// offsets between the hydrogens of a rotating methyl group, are captured
    /// whatever the angles' spread.
    #[default]
    Complex,
    /// Centered on the circular mean and wrapped to (-P/2, P/2], then treated
    /// as a real coordinate. Adequate for coordinates confined to one narrow
    /// well, but a correlation between coordinates that sweep the circle looks
    /// like a weak sawtooth and is largely missed.
    Centered,
}

/// Settings specific to the quasi-harmonic estimators.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct QuasiHarmonicOptions {
    /// How BQH estimates each coordinate's first-order entropy.
    pub marginals: MarginalEstimator,
    /// How periodic coordinates enter the correlation matrix.
    pub periodic: PeriodicCorrelation,
}

/// QH and BQH totals with the terms they are built from, all in units of k_B.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuasiHarmonicEntropy {
    /// `first_order_gaussian + correlation`. Without periodic coordinates this
    /// is 1/2 ln((2 pi e)^n |sigma|).
    pub quasi_harmonic: f64,
    /// `first_order_boltzmann + correlation`.
    pub boltzmann_quasi_harmonic: f64,
    /// Sum over coordinates of 1/2 ln(2 pi e sigma_ii).
    pub first_order_gaussian: f64,
    /// Sum over coordinates of -integral p ln p, from the marginal estimator.
    pub first_order_boltzmann: f64,
    /// 1/2 ln|C|, never positive.
    pub correlation: f64,
}

/// QH and BQH entropy of the coordinates in `one_d_data`.
///
/// `options.metrics` marks the periodic coordinates. The correlation matrix
/// must be nonsingular, which needs more frames than coordinates and no
/// coordinate that is a linear combination of others.
pub fn estimate_quasi_harmonic_entropy(
    one_d_data: &[Vec<f64>],
    frames_end: usize,
    quasi_harmonic: &QuasiHarmonicOptions,
    options: &EntropyOptions,
) -> Result<QuasiHarmonicEntropy, String> {
    let (one_d_data, metrics) = resolve_input(one_d_data, frames_end, options)?;
    let dimensions = one_d_data.len();
    if frames_end <= dimensions {
        return Err(format!(
            "the correlation matrix of {dimensions} coordinates needs more than {dimensions} frames, got {frames_end}"
        ));
    }

    let centered: Vec<Vec<f64>> = one_d_data
        .par_iter()
        .zip(metrics.par_iter())
        .map(|(values, metric)| center(values, *metric))
        .collect();
    let first_order_gaussian = centered
        .iter()
        .map(|values| {
            let variance = mean_square(values);
            0.5 * (2.0 * std::f64::consts::PI * std::f64::consts::E * variance).ln()
        })
        .sum();

    let deviations: Vec<Deviations> = centered
        .into_par_iter()
        .zip(one_d_data.par_iter().zip(metrics.par_iter()))
        .map(
            |(centered, (values, metric))| match (metric, quasi_harmonic.periodic) {
                (CoordinateMetric::Periodic { period }, PeriodicCorrelation::Complex) => {
                    complex_deviations(values, *period)
                }
                _ => Deviations {
                    real: centered,
                    imaginary: None,
                },
            },
        )
        .map(Deviations::standardized)
        .collect();
    let correlation = 0.5 * log_determinant(&correlation_matrix(&deviations))?;

    let first_order_boltzmann = match quasi_harmonic.marginals {
        MarginalEstimator::NearestNeighbor => coordinate_entropy_impl(
            &one_d_data,
            &metrics,
            options.constant,
            options.theiler_window,
        )?
        .iter()
        .sum(),
        MarginalEstimator::Histogram { bins } => {
            if bins == 0 {
                return Err("histogram needs at least one bin".to_string());
            }
            one_d_data
                .iter()
                .zip(&metrics)
                .map(|(values, metric)| histogram_entropy(values, *metric, bins))
                .sum()
        }
    };

    Ok(QuasiHarmonicEntropy {
        quasi_harmonic: first_order_gaussian + correlation,
        boltzmann_quasi_harmonic: first_order_boltzmann + correlation,
        first_order_gaussian,
        first_order_boltzmann,
        correlation,
    })
}

fn mean_square(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>() / values.len() as f64
}

/// Deviations of one coordinate from its mean, complex for a periodic
/// coordinate under [`PeriodicCorrelation::Complex`].
struct Deviations {
    real: Vec<f64>,
    imaginary: Option<Vec<f64>>,
}

impl Deviations {
    /// Scaled to unit mean |deviation|^2, so the correlation matrix has a unit
    /// diagonal.
    fn standardized(mut self) -> Self {
        let variance = mean_square(&self.real) + self.imaginary.as_deref().map_or(0.0, mean_square);
        let scale = variance.sqrt().recip();
        self.real.iter_mut().for_each(|v| *v *= scale);
        if let Some(imaginary) = self.imaginary.as_mut() {
            imaginary.iter_mut().for_each(|v| *v *= scale);
        }
        self
    }
}

/// z - <z> for z = e^(i 2 pi q / P).
fn complex_deviations(values: &[f64], period: f64) -> Deviations {
    let to_radians = std::f64::consts::TAU / period;
    let (real, imaginary): (Vec<f64>, Vec<f64>) = values
        .iter()
        .map(|v| (v * to_radians).cos())
        .zip(values.iter().map(|v| (v * to_radians).sin()))
        .unzip();
    let center = |series: Vec<f64>| {
        let mean = series.iter().sum::<f64>() / series.len() as f64;
        series.into_iter().map(|v| v - mean).collect::<Vec<f64>>()
    };
    Deviations {
        real: center(real),
        imaginary: Some(center(imaginary)),
    }
}

/// Deviations from the mean. A periodic coordinate is first rotated so its
/// circular mean sits at 0 and wrapped to (-P/2, P/2].
fn center(values: &[f64], metric: CoordinateMetric) -> Vec<f64> {
    let shifted: Vec<f64> = match metric {
        CoordinateMetric::Linear => values.to_vec(),
        CoordinateMetric::Periodic { period } => {
            let to_radians = std::f64::consts::TAU / period;
            let (sin, cos) = values.iter().fold((0.0, 0.0), |(s, c), v| {
                let angle = v * to_radians;
                (s + angle.sin(), c + angle.cos())
            });
            let circular_mean = sin.atan2(cos) / to_radians;
            let half = 0.5 * period;
            values
                .iter()
                .map(|v| half - (half - (v - circular_mean)).rem_euclid(period))
                .collect()
        }
    };
    let mean = shifted.iter().sum::<f64>() / shifted.len() as f64;
    shifted.iter().map(|v| v - mean).collect()
}

/// A complex number as (real, imaginary).
type Complex = (f64, f64);

/// C_ij = mean over frames of a_i conj(a_j) for standardized deviations a.
/// Hermitian with a unit diagonal; real when no deviation is complex.
fn correlation_matrix(deviations: &[Deviations]) -> Vec<Vec<Complex>> {
    let dimensions = deviations.len();
    let frames = deviations[0].real.len() as f64;
    let dot = |left: &[f64], right: &[f64]| left.iter().zip(right).map(|(a, b)| a * b).sum::<f64>();
    let upper: Vec<Vec<Complex>> = (0..dimensions)
        .into_par_iter()
        .map(|i| {
            (i..dimensions)
                .map(|j| {
                    if i == j {
                        return (1.0, 0.0);
                    }
                    let (a, b) = (&deviations[i], &deviations[j]);
                    // (ar + i ai)(br - i bi) = ar br + ai bi + i (ai br - ar bi)
                    let mut real = dot(&a.real, &b.real);
                    let mut imaginary = 0.0;
                    if let Some(ai) = &a.imaginary {
                        imaginary += dot(ai, &b.real);
                    }
                    if let Some(bi) = &b.imaginary {
                        imaginary -= dot(&a.real, bi);
                    }
                    if let (Some(ai), Some(bi)) = (&a.imaginary, &b.imaginary) {
                        real += dot(ai, bi);
                    }
                    (real / frames, imaginary / frames)
                })
                .collect()
        })
        .collect();
    let mut matrix = vec![vec![(0.0, 0.0); dimensions]; dimensions];
    for (i, row) in upper.iter().enumerate() {
        for (offset, &(real, imaginary)) in row.iter().enumerate() {
            matrix[i][i + offset] = (real, imaginary);
            matrix[i + offset][i] = (real, -imaginary);
        }
    }
    matrix
}

/// ln|A| of a Hermitian positive-definite matrix, from its Cholesky factor
/// A = L L^H, whose diagonal is real and positive.
fn log_determinant(matrix: &[Vec<Complex>]) -> Result<f64, String> {
    let n = matrix.len();
    let mut lower = vec![vec![(0.0, 0.0); n]; n];
    let mut log_det = 0.0;
    for j in 0..n {
        let diagonal = matrix[j][j].0
            - lower[j][..j]
                .iter()
                .map(|(re, im)| re * re + im * im)
                .sum::<f64>();
        // A NaN or non-positive pivot means C is not positive definite.
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(format!(
                "the correlation matrix is singular: coordinate {j} is a linear combination of earlier coordinates"
            ));
        }
        let pivot = diagonal.sqrt();
        lower[j][j] = (pivot, 0.0);
        log_det += 2.0 * pivot.ln();
        for i in (j + 1)..n {
            // sum over k of L_ik conj(L_jk)
            let (mut re, mut im) = (0.0, 0.0);
            for (a, b) in lower[i][..j].iter().zip(&lower[j][..j]) {
                re += a.0 * b.0 + a.1 * b.1;
                im += a.1 * b.0 - a.0 * b.1;
            }
            lower[i][j] = ((matrix[i][j].0 - re) / pivot, (matrix[i][j].1 - im) / pivot);
        }
    }
    Ok(log_det)
}

/// -sum p_k ln(p_k / width) over `bins` equal bins.
fn histogram_entropy(values: &[f64], metric: CoordinateMetric, bins: usize) -> f64 {
    let (low, width) = match metric {
        CoordinateMetric::Periodic { period } => (0.0, period / bins as f64),
        CoordinateMetric::Linear => {
            let low = values.iter().copied().fold(f64::INFINITY, f64::min);
            let high = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            (low, (high - low) / bins as f64)
        }
    };
    let mut counts = vec![0usize; bins];
    for &value in values {
        let value = metric.canonicalize(value);
        // The maximum of a linear range lands exactly on the upper edge.
        let bin = (((value - low) / width) as usize).min(bins - 1);
        counts[bin] += 1;
    }
    let total = values.len() as f64;
    counts
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / total;
            -p * (p / width).ln()
        })
        .sum()
}
