# nn-entropy

[![CI](https://github.com/helmutcarter/nn-entropy/actions/workflows/ci.yml/badge.svg)](https://github.com/helmutcarter/nn-entropy/actions/workflows/ci.yml)

Estimate the configurational entropy of a molecular system, with a Rust library, CLI, and Python bindings.

This crate provides non-parametric entropy estimation using the nearest-neighbor method with mutual information expansion on internal coordinates (bond lengths, bond angles, and torsion angles) with built-in conversion from MD trajectory files. File-based calculations use linear metrics for bonds and angles and a 2π-periodic metric for torsions.

## Features
- Rust library API for MIE entropy up to fourth order, per-coordinate entropy, and mutual information estimates.
- CLI that reads `.parm7` + `.nc` and converts to internal coordinates and prints total configurational entropy.
- Python bindings via `pyo3` for in-memory arrays or direct file-based calculation.
- Entropy calculations parallelized with `rayon`.

## Requirements
- Rust 1.94 or newer.
- For Python bindings: a Python environment with build tooling for `pyo3` (see below).

## Build

```bash
cargo build --release
```

## CLI usage

```bash
cargo run --release <path_to_parm7> <path_to_nc> [--torsions-only] [--no-periodic] [--exact-constant] [--start N] [--stop N] [--stride N] [--mie-order 1|2|3|4] [--theiler-window W] [--tree-bonds] [--no-phase-angles] [--jacobian] [--legacy] [--estimator nn|qh|bqh] [--histogram-bins N] [--periodic-correlation complex|centered]
```

Example:

```bash
cargo run --release <path_to_parm7> <path_to_nc>
```

Notes:
- `--stop` limits the number of frames read to N.
- `--start` skips the first N frames.
- `--stride` retains every Nth frame after `--start` (default: 1). Choose a stride based on a separate correlation-time analysis; the crate does not estimate an effective sample size automatically.
- `--mie-order` selects the expansion order. The default is 2, and this should be good for most applications.
- `--no-periodic` forces every coordinate to use an ordinary linear metric. This produces a significant speedup, but the torsion entropy contribution loses some correctness.
- `--theiler-window W` excludes temporally adjacent frames from the nearest-neighbor search: frame i's neighbor must satisfy |i - j| > W, counted in frames after `--start` and `--stride`. The default, 0, is the plain estimator. See [Correlated frames](#correlated-frames-and-the-theiler-window).
- `--exact-constant` uses the exact finite-sample term `psi(N) - psi(1) = H_(N-1)`. The default is the common approximation `ln(N) + EulerGamma`.
- `--tree-bonds` uses the bonds of the BAT spanning tree, including bonds to hydrogen and excluding ring closures, in place of every heavy-atom bond. See [Coordinate sets](#coordinate-sets).
- `--no-phase-angles` uses the raw torsions. By default each torsion that shares its central bond with an earlier torsion is replaced by its offset from that torsion. See [Coordinate sets](#coordinate-sets).
- `--jacobian` also prints `Mean ln J = ...`, the mean of the log BAT Jacobian over the frames used, on the line before the entropy. It requires `--tree-bonds` and cannot be combined with `--torsions-only`.
- `--legacy` reproduces the original Python implementation bit for bit. See [Legacy mode](#legacy-mode).
- `--estimator` chooses the estimator: `nn` (the default) is the nearest-neighbor MIE, `qh` is quasi-harmonic, and `bqh` is Boltzmann quasi-harmonic. See [Quasi-harmonic estimators](#quasi-harmonic-estimators). `--mie-order` applies only to `nn`.
- `--histogram-bins N` makes `bqh` estimate each coordinate's first-order entropy from an N-bin histogram instead of nearest neighbors.
- `--periodic-correlation complex|centered` sets how torsions enter the quasi-harmonic correlation term; `complex`, the default, is Harpole and Sharp's `e^{iφ}`. See [Quasi-harmonic estimators](#quasi-harmonic-estimators).

## Coordinate sets

By default, file-based calculations use every bond between two heavy atoms (ring closures included), plus the BAT tree's angles and torsions, with phase angles in place of redundant torsions. Bonds to hydrogen are left out because they are usually constrained.

**Tree bonds** (`--tree-bonds`, `tree_bonds=True`, `BondSet::Tree`) use the bonds of the BAT spanning tree: the two root bonds and the bond that places each later atom, with hydrogens included and ring closures left out. Together with the angles and torsions, this gives exactly 3N - 6 non-redundant coordinates per molecule (one bond for a diatomic, including one with hydrogen). This is the coordinate set the BAT Jacobian is defined for.

**Phase angles** (on by default; `--no-phase-angles`, `phase_angles=False` or `CoordinateOptions::phase_angles` turn them off) group the torsions by central bond. The first torsion about each bond, in coordinate order, keeps its value. Every later torsion about that bond is replaced by its offset from the first one, wrapped to (-π, π]. The offset is still a 2π-periodic coordinate. Because the map is a shear with unit Jacobian, the joint entropy of the full coordinate set is unchanged. A truncated MIE is not invariant, though. The raw torsions of atoms that rotate together, such as the three hydrogens of a methyl group, are strongly correlated, and a second-order expansion misses part of that correlation. Their offsets are stiff and nearly independent of the reference torsion, so less correlation is left beyond second order. Every estimator in the crate does better with them; see [Validation](#validation). Phase angles work with either bond set and with `--torsions-only`.

**The BAT Jacobian** (`--jacobian`, `System.mean_log_jacobian`, `InternalCoordinates::log_jacobian`) is ln J = Σ 2 ln b + Σ ln sin θ over the tree's bonds and angles. Adding its mean to the BAT entropy gives the entropy of the internal motion measured in Cartesian coordinates. Bond lengths are in the trajectory's units (Å for Amber), so ln J depends on that choice, just as the bond-length entropies do. Phase angles do not change it.

## Quasi-harmonic estimators

Both quasi-harmonic estimators describe correlations with the correlation matrix C of the coordinates, `C_ij = σ_ij / √(σ_ii σ_jj)`, which contributes `½ ln|C| ≤ 0`. They differ in the first-order term:

- **QH** (`--estimator qh`, Karplus and Kushick 1981) fits a Gaussian to each coordinate: `S = Σ ½ ln(2πe σ_ii) + ½ ln|C|`, which is `½ ln((2πe)^n |σ|)` when no coordinate is periodic.
- **BQH** (`--estimator bqh`, Harpole and Sharp, *J. Phys. Chem. B* 2011) keeps the correlation term but takes each coordinate's first-order entropy, `-∫ p ln p`, from its actual distribution, so multimodal coordinates are not treated as one wide well.

The BQH first-order term uses the one-dimensional nearest-neighbor estimate by default, the same first-order term as the MIE, so BQH and the NN MIE differ only in how they treat correlations; `--exact-constant` and `--theiler-window` apply to it. `--histogram-bins N` uses N equal bins instead: linear coordinates are binned over their observed range, periodic ones over a full period. The histogram gives the differential entropy `-Σ p_k ln(p_k / Δ)`, Harpole and Sharp's eq 9, comparable with the other estimators. Some later BQH work reports the discrete `-Σ p_k ln p_k`, which is this value minus `ln Δ` per coordinate, a difference that cancels in entropy changes computed with the same bins. Harpole and Sharp choose bins by width (3-5° for angles, 0.1 Å for lengths) rather than by count. All values are differential entropies in k_B, like the nearest-neighbor results.

The covariance is the population covariance over the selected frames, so there must be more frames than coordinates, and no coordinate can be a linear combination of others.

**Torsions.** By default (`--periodic-correlation complex`), each periodic coordinate enters the correlation matrix as `z = e^{iφ}` (for period P, `e^{2πiq/P}`), following Harpole and Sharp, which makes C Hermitian; `C_ij = ⟨δa_i δa_j*⟩ / √(⟨|δa_i|²⟩⟨|δa_j|²⟩)` with `δa = a - ⟨a⟩`, and linear coordinates stay real. This captures correlations on the circle whatever the angles' spread, such as the fixed 120° offsets between the three hydrogen torsions of a rotating methyl group. The first-order Gaussian term of QH uses each torsion's variance after rotating its circular mean to 0 and wrapping to (-π, π].

`--periodic-correlation centered` puts the centered and wrapped torsions into C as real coordinates instead. That is adequate for torsions confined to one narrow well, where both treatments agree, but a correlation between torsions that sweep the circle looks like a weak sawtooth after wrapping and is largely missed.

See [Validation](#validation) for how the estimators compare with thermodynamic reference values.

## Validation

The estimators have been checked against thermodynamic reference entropies: on the butane end states, BQH with phase angles agrees with thermodynamic integration within the reference's uncertainty. The [validation page](https://helmutcarter.github.io/nn-entropy/nn_entropy/validation/index.html) of the API documentation has the results for every estimator and an example program that reproduces them.

## Legacy mode

`--legacy` gives the same total, to the last bit, as the original Python script `NN_entropy_calc.py`. Use it to reproduce historical results or to compare against them. It is not more accurate than the default; it undoes the corrections made since. It changes four things:

- Every coordinate uses a linear metric, torsions included. This is the only change larger than rounding.
- The estimator follows the script's arithmetic: order-2 MIE, gamma truncated to `0.57721566490153`, constants `ln(2N) + gamma` and `ln(πN) + gamma`, and per-coordinate and per-pair entropies added with Python's built-in `sum`, which uses compensated summation from Python 3.12 on.
- Angles use the script's `acos` formula (`Geometry::Legacy`), which can return NaN for a straight angle.
- Torsions take |b2| as `pow(x, 0.5)`, as NumPy evaluates `x ** 0.5`, which is off by one ulp from `sqrt(x)` for about 1 in 1000 inputs.

It uses the raw torsions, as the script did, and cannot be combined with `--exact-constant`, `--theiler-window`, `--mie-order` other than 2, `--estimator` other than `nn`, `--tree-bonds`, `--phase-angles` or `--jacobian`, because the script had none of them. `--torsions-only`, `--start`, `--stop` and `--stride` only choose the input and are allowed. The script's `N` argument corresponds to `--stop N`.

Atom selection is this crate's, which handles layouts the script did not (see [Interpretation and limitations](#interpretation-and-limitations)). The script took the first atoms of the trajectory and recognized water only by the types `OW`, `HW` and `EP`. On systems it handled correctly, both select the same atoms. Bit-for-bit agreement was checked with Python 3.12, NumPy 1.26, SciPy 1.13 and MDAnalysis 2.8 on seven trajectories: the test fixture, five cpptraj sample systems (up to 821 coordinates, including a solvated one), and two 10,000-frame butane end states. A different Python or NumPy version can change the script's result in the last bits. From Rust, the mode is `bat_library::Geometry::Legacy` for the coordinates and `legacy::legacy_entropy` for the estimator.

## Rust library usage

```rust
use nn_entropy::{CoordinateMetric, EntropyOptions, FiniteSampleConstant, calculate_entropy};

// one_d_data is Vec<Vec<f64>> with shape [n_coords][n_frames]
let entropy = calculate_entropy(&one_d_data, frames_end, 2, &EntropyOptions::default())?;

// Periodic torsions and the exact finite-sample constant
let options = EntropyOptions {
    metrics: Some(&metrics), // one CoordinateMetric per coordinate
    constant: FiniteSampleConstant::Exact,
    ..Default::default()
};
let entropy = calculate_entropy(&one_d_data, frames_end, 2, &options)?;
```

`EntropyOptions` has three fields:
- `metrics`: per-coordinate `CoordinateMetric::Linear` or `CoordinateMetric::Periodic { period }`. The default, `None`, treats every coordinate as linear.
- `constant`: `FiniteSampleConstant::Asymptotic` (the default) or `FiniteSampleConstant::Exact`.
- `theiler_window`: frames j with |i - j| <= `theiler_window` are never frame i's nearest neighbor. The default, 0, excludes only frame i itself.

Every estimator takes the same options:
- `calculate_entropy(data, frames_end, mie_order, &options)` returns the total entropy at MIE order 1, 2, 3 or 4.
- `estimate_coordinate_entropy(data, frames_end, mie_order, &options)` returns each coordinate's share of the entropy at MIE order 1 or 2, so the values sum to `calculate_entropy` at the same order and with the same options. Order 1 gives the marginal entropies. Order 2 also splits each pairwise mutual information term evenly between its two coordinates. Orders 3 and 4 are rejected.
- `estimate_coordinate_mutual_information(data, frames_end, &options)` returns the mutual information of each coordinate pair.
- `calc_joint_nn([&x, &y, ...], metrics)` is the underlying nearest-neighbor sum, `sum over frames of ln(nearest-neighbor distance)`, in 1 to 4 dimensions. `calc_joint_nn_with_window(coordinates, metrics, w)` is the same sum with a Theiler window.

To build coordinates from files, use `bat_library::InternalCoordinates`:

```rust
use nn_entropy::bat_library::{BondSet, CoordinateOptions, InternalCoordinates};

let options = CoordinateOptions { bonds: BondSet::Tree, ..Default::default() }; // phase angles by default
let mut internal = InternalCoordinates::with_options(&top, options)?; // ::new(&top) for the defaults
internal.calculate_internal_coords(&traj, usize::MAX, false)?;
let metrics = internal.coordinate_metrics();   // one per coordinate
let atoms = internal.coordinate_atoms();       // topology atom indices of each coordinate
let references = internal.phase_references();  // Some(i) for a phase angle of torsion i
let log_jacobian = internal.log_jacobian()?;   // per frame; tree bonds only
```

`internal.int_coords` holds the values with shape `[n_frames][n_coords]`.

The quasi-harmonic estimators take the same data and options:

```rust
use nn_entropy::quasi_harmonic::{
    MarginalEstimator, PeriodicCorrelation, QuasiHarmonicOptions, estimate_quasi_harmonic_entropy,
};

let quasi_harmonic = QuasiHarmonicOptions::default(); // NN marginals, complex torsions
// or QuasiHarmonicOptions { marginals: MarginalEstimator::Histogram { bins: 100 }, periodic: PeriodicCorrelation::Centered }
let entropy = estimate_quasi_harmonic_entropy(&one_d_data, frames_end, &quasi_harmonic, &options)?;
println!("QH {} BQH {}", entropy.quasi_harmonic, entropy.boltzmann_quasi_harmonic);
```

The result also carries the terms the totals are built from: `first_order_gaussian`, `first_order_boltzmann` and `correlation`.

## Python bindings

The crate exposes a `nn_entropy` Python module (built from `src/pyo3_api.rs`) with:
- `load_system(top_path, traj_path, start=None, stop=None, torsions_only=None, stride=None, tree_bonds=False, phase_angles=True)`. See [Coordinate sets](#coordinate-sets).
- `estimate_entropy(data_or_system, mie_order=None, periods=None, exact_constant=False, theiler_window=0)`. `mie_order` defaults to 2; `mie_order=1` gives the sum of the marginal entropies.
- `estimate_coordinate_entropy(data_or_system, mie_order=None, periods=None, exact_constant=False, theiler_window=0)`. Returns each coordinate's share of the entropy; `mie_order` defaults to 2 and must be 1 or 2. The values sum to `estimate_entropy` at the same order.
- `estimate_coordinate_mutual_information(data_or_system, periods=None, exact_constant=False, theiler_window=0)`. Returns the mutual information of each coordinate pair.
- `estimate_quasi_harmonic_entropy(data_or_system, periods=None, histogram_bins=None, periodic_correlation="complex", exact_constant=False, theiler_window=0)`. Returns a dict with `quasi_harmonic`, `boltzmann_quasi_harmonic`, `first_order_gaussian`, `first_order_boltzmann` and `correlation`. See [Quasi-harmonic estimators](#quasi-harmonic-estimators).
- `System` has the same four estimators as methods, each taking `exact_constant` and `theiler_window` (and `mie_order` for the two MIE entropy estimators, `histogram_bins` and `periodic_correlation` for the quasi-harmonic one). Its `mean_log_jacobian` property is the mean log BAT Jacobian over the selected frames. It raises `ValueError` unless the system was loaded with `tree_bonds=True` and without `torsions_only`.
- To start from an Amber topology and trajectory, call `load_system(...)` once and use the `System` methods or pass the `System` to any estimator. The coordinate metrics come from the generated BAT coordinates automatically, and the trajectory is only read once no matter how many quantities you compute.

`exact_constant=True` uses the exact finite-sample term, like the CLI's `--exact-constant`, and `theiler_window` matches the CLI's `--theiler-window`.

For a raw array, `periods` is a sequence with one entry per coordinate. Use `None` for a linear dimension and a positive period for a periodic dimension. Omitting `periods` preserves the linear behavior of the original raw-array API.

Preferred usage:

```python
import nn_entropy

system = nn_entropy.load_system("system.parm7", "trajectory.nc")
entropy = nn_entropy.estimate_entropy(system, mie_order=2)

coordinate_entropy = system.estimate_coordinate_entropy()  # order 2, sums to `entropy`
marginal_entropy = system.estimate_coordinate_entropy(mie_order=1)
```

A typical build workflow uses `maturin`:

```bash
maturin develop --release
```

Threading note:
- The Python wrappers release the GIL while loading trajectories and running entropy calculations, so multiple Python threads can call into `nn_entropy` concurrently.
- Internal parallelism is handled in Rust with Rayon. Set `RAYON_NUM_THREADS=N` if you want to cap or tune the Rust worker pool.

## Interpretation and limitations

- File-based calculations select the solute automatically. Water is excluded when its atom type is an Amber water type (`OW`, `HW`, `EP`, `EPW`) or its residue is a water name (`WAT`, `HOH`, `SOL`, `TIP3`, `SPCE`, `OPC`, and similar), so four-point models and non-Amber water types are handled. Massless virtual sites are excluded because they carry no degrees of freedom. Monatomic species such as counter-ions contribute nothing, since they have no internal coordinates, and a diatomic contributes only its bond length (with the default bond set, only when both atoms are heavy). Solute atoms may appear anywhere in the topology.
- The trajectory must have exactly as many atoms as the topology. A stripped topology paired with a solvated trajectory is rejected rather than read on the assumption that atoms line up.
- Results are differential internal-coordinate entropies in natural-log units (equivalently, units of `k_B`), not kcal mol⁻¹ K⁻¹.
- Following the historical Python implementation, the default nearest-neighbor constant uses the large-sample approximation `ln(N) + EulerGamma` in place of the exact finite-sample term `psi(N) - psi(1) = H_(N-1)`. Their per-entropy difference is approximately `1/(2N)` and is negligible for the intended 50,000-frame calculations, but it can accumulate across MIE terms or become important for very small samples. The convention is a property of the estimator, not of the total alone: pairwise mutual information shifts by exactly `H_(N-1) - (ln(N) + EulerGamma)` per pair when it changes, so per-coordinate results must be computed under the same convention as the total they are compared against. The Python implementation truncated EulerGamma to `0.57721566490153`; this crate uses the full double-precision value, so results can differ from it by roughly 1e-14. `--legacy` uses the truncated value.
- The default second-order mutual-information expansion is a truncation and can omit higher-order correlations.
- The reported value does not include a BAT-to-Cartesian Jacobian correction, momentum entropy, or a standard-state term, and must not be described as an absolute thermodynamic entropy. `--jacobian` reports the Jacobian term separately, for the tree bond set, but does not add it to the total.
- Exact duplicate samples use the nearest distinct coordinate point. Many ties usually indicate inadequate coordinate precision or sampling and should be investigated.
- MD frames are generally correlated. See below.
- File-based results differ from the historical Python implementation because torsions are now periodic and phase angles replace redundant torsions. Its asymptotic finite-sample convention is kept as the default. `--legacy` reproduces the historical results exactly.

## Correlated frames and the Theiler window

The Kozachenko-Leonenko estimator assumes independent samples. When frames are saved faster than a coordinate decorrelates, a frame's nearest neighbor is usually the frame before or after it. That distance measures how far the system moves in one output interval, not the density of states, so entropies come out low. The error grows with the number of dimensions, so it does not cancel in the MIE: oversampled independent coordinates show spurious positive mutual information.

A Theiler window of `W` frames keeps every frame as a query point but takes its nearest neighbor only from frames more than `W` away in time. Each frame's finite-sample term uses the number of frames it can actually compete against, so the window does not bias the estimate for independent samples. Unlike `--stride`, it discards no data. Frames are counted after `--start` and `--stride`, and the series must have at least `2W + 2` frames.

What the validation shows (`tests/theiler.rs`):
- On an AR(1) series with lag-one correlation 0.99, the mutual information of two independent coordinates is +0.044 at `W = 0` and about 0 at `W = 5`.
- On mildly correlated series (lag-one correlation 0.95), any `W` of a few frames or more removes the bias and leaves the estimate on a plateau. RMSE is lower than with any stride.
- On strongly oversampled series the window over-corrects as `W` grows. The frames that remain eligible come in clusters, one per visit to a region, which lengthens nearest-neighbor distances. Mutual information then turns negative, about -0.04 at `W >= 50` in the example above. A window is not a substitute for adequate sampling.
- On 100 ns production runs of small solvated ligands saved every 2 ps (50,000 frames), the order-2 total changes by at most about 0.4 k_B across `W = 0` to `250` (three ligands), while striding moves it by 10-20 k_B. At that spacing the frames are already effectively decorrelated, so the plain estimator is fine there.

Use a window when output is dense compared with the coordinates' correlation times. Choose `W` from a sweep: plot the entropy against `W` and take a value where the curve has stopped moving for the fast coordinates. Prefer a window to a stride when the goal is only to remove temporal neighbors.

## Documentation

The Rust API documentation, including a guide to choosing an estimator and the validation results, is published at <https://helmutcarter.github.io/nn-entropy/> for the current master. To build it locally:

```bash
cargo doc --no-deps --open
```

## Tests

```bash
cargo test --release
```

The suite is split by concern:

| File | Covers |
|---|---|
| `tests/common/` | Shared helpers: a self-contained RNG, closed-form entropies, an O(N^2) nearest-neighbor reference, and writers for Amber topologies and NetCDF trajectories (`amber.rs`). Calls none of the code it checks. |
| `tests/analytic.rs` | Closed-form ground truth: uniform, normal, exponential, circular uniform, bivariate-normal joint entropy and mutual information, plus a convergence check. |
| `tests/reference.rs` | Brute-force cross-checks of the kd-tree search for 1-4 dimensions and every metric combination, and of the MIE expansion against an independent inclusion-exclusion evaluation. |
| `tests/invariance.rs` | Translation, reflection, relabeling, frame reordering, scaling equivariance, periodic rotation, and floating-point behavior. |
| `tests/geometry.rs` | Bonds, angles and torsions from constructed configurations; the torsion sign convention; the end-to-end BAT path on the fixture. |
| `tests/topology.rs` | Atom selection and file reading across layouts the fixture does not cover: solvent before or between solute molecules, counter-ions, four-point and non-Amber water, and NetCDF corner cases. |
| `tests/validation.rs` | The input-validation contract: which inputs are refused and what each error says. |
| `tests/coordinates.rs` | Coordinate sets: 3N - 6 tree bonds forming a spanning tree, phase angles as wrapped torsion offsets that are stiff for a rigid rotor, and the BAT Jacobian against the constructed geometry. |
| `tests/quasi_harmonic.rs` | QH against the Gaussian entropy of independently computed sample covariances and the bivariate-normal closed form, BQH correcting QH for a bimodal coordinate, histograms against hand-counted bins, the complex torsion correlations against hand-computed Hermitian coefficients and determinants, periodic centering, validation, and the CLI flags. |
| `tests/legacy.rs` | `--legacy` against totals produced by the original Python script on the fixture, compared bit for bit, and the legacy angle and torsion formulas. |
| `tests/theiler.rs` | Theiler window: bitwise equality with the plain estimator at `W = 0`, brute-force agreement for every window and dimension, the closed-form constant, validation, and AR(1) bias behavior. |
| `tests/testing.rs`, `tests/entropy.rs`, `tests/cli.rs` | Unit tests, fixture regressions, and CLI behavior. |

Expected values in `analytic.rs`, `reference.rs`, `invariance.rs`, `geometry.rs`,
`coordinates.rs` and `theiler.rs` are derived from theory or from an independent
reference, never recorded from a previous run. The pinned totals in `entropy.rs`
are regression values.

The Rust suite does not link the Python extension, so the bindings have a
separate smoke test that runs against a built wheel:

```bash
maturin build --out dist
pip install --no-index --find-links dist nn_entropy
python tests/python_smoke.py
```

CI (`.github/workflows/ci.yml`) runs `cargo fmt --check`, `cargo clippy -D warnings`
and `cargo test` on every push and pull request, plus the wheel build and smoke
test on Python 3.10 and 3.12. A separate job checks that the crate still builds
on the minimum supported Rust version, 1.94 (`rust-version` in `Cargo.toml`).

## Project layout
- `src/lib.rs`: the nearest-neighbor estimators, the MIE, and the nearest-neighbor search; also the crate documentation.
- `src/quasi_harmonic.rs`: the QH and BQH estimators.
- `src/legacy.rs`: the reproduction of the original Python implementation.
- `src/validation.rs`: documentation only, the validation results.
- `src/bat_library/`: `.parm7` parser, NetCDF trajectory reader, solute selection, and BAT coordinate construction (`geometry.rs` holds the bond, angle and torsion functions).
- `src/main.rs`: CLI entry point.
- `src/pyo3_api.rs`: Python bindings.
- `tests/`: the test suite (see [Tests](#tests)).
- `examples/butane_validation.rs`: reproduces the validation results from the butane trajectories.

## Further Reading
- [Grid inhomogeneous solvation theory: Hydration structure and thermodynamics of the miniature receptor cucurbit[7]uril](https://pmc.ncbi.nlm.nih.gov/articles/PMC3416872/) - Uses nearest neighbor method to calculate first-order estimate of solvent entropy
- [Extraction of configurational entropy from molecular simulations via an expansion approximation](https://pubmed.ncbi.nlm.nih.gov/17640119/) - Uses mutual information expansion to increase accuracy of entropy estimation of highly correlated systems like molecules
- [Sample Estimate of the Entropy of a Random Vector](https://dmitripavlov.org/scans/kozachenko-leonenko.pdf) - First introduction of nearest neighbors entropy estimation


## License
© 2026 Helmut Carter, Kurtzman Lab. All rights reserved.
