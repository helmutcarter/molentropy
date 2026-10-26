//! How an estimate changes with the number of frames.
//!
//! Every estimator in the crate is biased at a finite number of frames: the
//! nearest-neighbor terms by an amount that grows with their dimension, and
//! Gaussian mutual information terms by about pq / 2n for p and q correlated
//! columns. [`check_convergence`] repeats an estimate on every 4th, every 2nd
//! and every frame. Thinning keeps the span of the trajectory, and so the
//! conformations it sampled, and changes only how many frames describe them,
//! so a change between the three estimates is finite-sample bias, not drift in
//! the sampling. It cannot detect states the trajectory never visited.
//!
//! The extrapolation assumes the bias falls as 1/n, as it does for Gaussian
//! terms and roughly for one- and two-dimensional nearest-neighbor terms. A
//! bias that falls more slowly, such as that of a joint nearest-neighbor
//! estimate in many dimensions, is only partly removed.

/// The strides of [`check_convergence`], from fewest frames to all of them.
pub const CONVERGENCE_STRIDES: [usize; 3] = [4, 2, 1];

/// One estimate of [`check_convergence`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvergencePoint {
    /// Every `stride`-th frame was used.
    pub stride: usize,
    /// Frames used.
    pub frames: usize,
    /// The estimate from them.
    pub entropy: f64,
}

/// The result of [`check_convergence`].
#[derive(Debug, Clone, PartialEq)]
pub struct Convergence {
    /// Estimates from every 4th, every 2nd and every frame.
    pub points: [ConvergencePoint; 3],
    /// The estimate from all frames minus that from every 2nd frame.
    pub change: f64,
    /// The estimate extrapolated to infinitely many frames; see
    /// [`extrapolate`].
    pub extrapolated: f64,
}

impl Convergence {
    /// Whether the estimate moves the same way from n/4 to n/2 frames as from
    /// n/2 to n. If not, sampling noise is larger than the bias, or the bias
    /// does not fall steadily, and the extrapolation means little.
    pub fn steady(&self) -> bool {
        let [quarter, half, full] = self.points.map(|point| point.entropy);
        (half - quarter) * (full - half) >= 0.0
    }
}

/// Every `stride`-th frame of the first `frames_end` frames of each coordinate.
pub fn subsample(one_d_data: &[Vec<f64>], frames_end: usize, stride: usize) -> Vec<Vec<f64>> {
    one_d_data
        .iter()
        .map(|column| {
            column[..frames_end]
                .iter()
                .step_by(stride)
                .copied()
                .collect()
        })
        .collect()
}

/// The intercept of the least-squares line through the estimates `entropies`
/// against 1 / `frames`: their limit for infinitely many frames if the bias
/// falls as 1/n. A fit through all points rather than an extrapolation from
/// successive differences, which noise in the estimates makes unstable.
pub fn extrapolate(frames: &[usize], entropies: &[f64]) -> f64 {
    let x: Vec<f64> = frames.iter().map(|&n| (n as f64).recip()).collect();
    let count = x.len() as f64;
    let x_mean = x.iter().sum::<f64>() / count;
    let y_mean = entropies.iter().sum::<f64>() / count;
    let covariance: f64 = x
        .iter()
        .zip(entropies)
        .map(|(x, y)| (x - x_mean) * (y - y_mean))
        .sum();
    let variance: f64 = x.iter().map(|x| (x - x_mean).powi(2)).sum();
    y_mean - covariance / variance * x_mean
}

/// Run `estimate` on every 4th, every 2nd and every frame of the first
/// `frames_end` frames. `estimate` receives the thinned coordinates, their
/// number of frames, and the stride; a Theiler window counts frames of the
/// thinned data, so a window of W frames becomes ceil(W / stride).
pub fn check_convergence<F>(
    one_d_data: &[Vec<f64>],
    frames_end: usize,
    mut estimate: F,
) -> Result<Convergence, String>
where
    F: FnMut(&[Vec<f64>], usize, usize) -> Result<f64, String>,
{
    if one_d_data.is_empty() {
        return Err("no coordinates".to_string());
    }
    if let Some(short) = one_d_data.iter().find(|column| column.len() < frames_end) {
        return Err(format!(
            "{frames_end} frames requested but a coordinate has {}",
            short.len()
        ));
    }
    let mut points = Vec::with_capacity(CONVERGENCE_STRIDES.len());
    for stride in CONVERGENCE_STRIDES {
        let data = subsample(one_d_data, frames_end, stride);
        let frames = data[0].len();
        let entropy = estimate(&data, frames, stride)
            .map_err(|err| format!("{}: {err}", describe(stride)))?;
        points.push(ConvergencePoint {
            stride,
            frames,
            entropy,
        });
    }
    let points: [ConvergencePoint; 3] = points.try_into().expect("three strides");
    let entropies = points.map(|point| point.entropy);
    Ok(Convergence {
        points,
        change: entropies[2] - entropies[1],
        extrapolated: extrapolate(&points.map(|point| point.frames), &entropies),
    })
}

fn describe(stride: usize) -> String {
    match stride {
        1 => "all frames".to_string(),
        2 => "every 2nd frame".to_string(),
        n => format!("every {n}th frame"),
    }
}
