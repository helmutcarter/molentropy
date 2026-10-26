//! The estimator of the original Python implementation, `NN_entropy_calc.py`,
//! reproduced operation for operation so the results agree to the last bit.
//!
//! It differs from [`crate::calculate_entropy`] at order 2 in ways that are
//! each tiny but change the floating-point result:
//! - every coordinate uses a linear metric, torsions included;
//! - gamma is the truncated 0.57721566490153;
//! - the constants are ln(2N) + gamma and ln(pi N) + gamma, and each entropy
//!   is `sum / N + constant` rather than `sum * (1 / N) + constant`;
//! - the per-coordinate and per-pair entropies are summed in order with
//!   Python's built-in `sum`, which is compensated from Python 3.12 on.
//!
//! Only the linear metric changes the result by more than rounding; see
//! <https://helmutcarter.github.io/molentropy/legacy.html>.

use crate::{CoordinateMetric, EntropyOptions, binomial, calc_joint_nn, combination_from_rank};
use rayon::prelude::*;

/// Gamma as `NN_entropy_calc.py` hard-codes it.
#[allow(clippy::approx_constant)]
const PYTHON_EULER_GAMMA: f64 = 0.57721566490153;

/// Order-2 MIE entropy exactly as `NN_entropy_calc.py` computes it.
pub fn legacy_entropy(one_d_data: &[Vec<f64>], frames_end: usize) -> Result<f64, String> {
    let (one_d_data, _) = crate::resolve_input(one_d_data, frames_end, &EntropyOptions::default())?;
    let frames = frames_end as f64;
    let one_d_constant = ((2 * frames_end) as f64).ln() + PYTHON_EULER_GAMMA;
    let two_d_constant = (frames * std::f64::consts::PI).ln() + PYTHON_EULER_GAMMA;
    let linear = CoordinateMetric::Linear;

    let one_d_entropies = one_d_data
        .par_iter()
        .map(|coordinate| {
            calc_joint_nn([coordinate], [linear]).map(|total| total / frames + one_d_constant)
        })
        .collect::<Result<Vec<f64>, String>>()?;

    let degrees_freedom = one_d_data.len();
    let two_d_entropies = (0..binomial(degrees_freedom, 2))
        .into_par_iter()
        .map(|rank| {
            let [i, j] = combination_from_rank::<2>(rank, degrees_freedom);
            calc_joint_nn([&one_d_data[i], &one_d_data[j]], [linear, linear])
                .map(|total| total * 2.0 / frames + two_d_constant)
        })
        .collect::<Result<Vec<f64>, String>>()?;

    Ok(
        python_sum(&two_d_entropies)
            - (degrees_freedom as f64 - 2.0) * python_sum(&one_d_entropies),
    )
}

/// Python's built-in `sum` over floats from 3.12 on: Neumaier's compensated
/// summation, as CPython implements it.
fn python_sum(values: &[f64]) -> f64 {
    let mut total = 0.0_f64;
    let mut compensation = 0.0_f64;
    for &value in values {
        let next = total + value;
        if total.abs() >= value.abs() {
            compensation += (total - next) + value;
        } else {
            compensation += (value - next) + total;
        }
        total = next;
    }
    if compensation != 0.0 && compensation.is_finite() {
        total += compensation;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::python_sum;

    #[test]
    fn python_sum_is_compensated() {
        // Python 3.12: sum([0.1] * 10) == 1.0, where a plain loop gives
        // 0.9999999999999999.
        assert_eq!(python_sum(&[0.1; 10]), 1.0);
        assert_eq!(python_sum(&[1e100, 1.0, -1e100]), 1.0);
        assert_eq!(python_sum(&[]), 0.0);
    }
}
