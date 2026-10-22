# Legacy mode

`--legacy` gives the same total, to the last bit, as the original Python script `NN_entropy_calc.py`. Use it to reproduce historical results or to compare against them. It is not more accurate than the default; it undoes the corrections made since. It changes four things:

- Every coordinate uses a linear metric, torsions included. This is the only change larger than rounding.
- The estimator follows the script's arithmetic: order-2 MIE, gamma truncated to `0.57721566490153`, constants `ln(2N) + gamma` and `ln(πN) + gamma`, and per-coordinate and per-pair entropies added with Python's built-in `sum`, which uses compensated summation from Python 3.12 on.
- Angles use the script's `acos` formula (`Geometry::Legacy`), which can return NaN for a straight angle.
- Torsions take |b2| as `pow(x, 0.5)`, as NumPy evaluates `x ** 0.5`, which is off by one ulp from `sqrt(x)` for about 1 in 1000 inputs.

It uses the raw torsions, as the script did, and cannot be combined with `--exact-constant`, `--theiler-window`, `--mie-order` other than 2, `--estimator` other than `nn`, `--tree-bonds`, `--phase-angles` or `--jacobian`, because the script had none of them. `--torsions-only`, `--start`, `--stop` and `--stride` only choose the input and are allowed. The script's `N` argument corresponds to `--stop N`.

Atom selection is this crate's, which handles layouts the script did not (see [Interpretation and limitations](limitations.md)). The script took the first atoms of the trajectory and recognized water only by the types `OW`, `HW` and `EP`. On systems it handled correctly, both select the same atoms. Bit-for-bit agreement was checked with Python 3.12, NumPy 1.26, SciPy 1.13 and MDAnalysis 2.8 on seven trajectories: the test fixture, five cpptraj sample systems (up to 821 coordinates, including a solvated one), and two 10,000-frame butane end states. A different Python or NumPy version can change the script's result in the last bits. From Rust, the mode is `bat_library::Geometry::Legacy` for the coordinates and `legacy::legacy_entropy` for the estimator.
