# Command line

```bash
cargo run --release <path_to_parm7> <path_to_nc> [--torsions-only] [--no-periodic] [--exact-constant] [--start N] [--stop N] [--stride N] [--mie-order 1|2|3|4] [--theiler-window W] [--tree-bonds] [--no-phase-angles] [--jacobian] [--legacy] [--estimator nn|qh|bqh] [--histogram-bins N] [--periodic-correlation complex|centered]
```

Example:

```bash
cargo run --release <path_to_parm7> <path_to_nc>
```

## Options

- For most molecules `--estimator bqh` is the most accurate choice; see [Validation](validation.md). The default estimator is the NN MIE.
- `--torsions-only` drops the bonds and angles and estimates the entropy of the torsions alone.
- `--stop` limits the number of frames read to N.
- `--start` skips the first N frames.
- `--stride` retains every Nth frame after `--start` (default: 1). Choose a stride based on a separate correlation-time analysis; the crate does not estimate an effective sample size automatically.
- `--mie-order` selects the expansion order of the NN MIE. The default is 2.
- `--no-periodic` forces every coordinate to use an ordinary linear metric. This produces a significant speedup, but the torsion entropy contribution loses some correctness.
- `--theiler-window W` excludes temporally adjacent frames from the nearest-neighbor search: frame i's neighbor must satisfy |i - j| > W, counted in frames after `--start` and `--stride`. The default, 0, is the plain estimator. See [Correlated frames](correlated-frames.md).
- `--exact-constant` uses the exact finite-sample term `psi(N) - psi(1) = H_(N-1)`. The default is the common approximation `ln(N) + EulerGamma`.
- `--tree-bonds` uses the bonds of the BAT spanning tree, including bonds to hydrogen and excluding ring closures, in place of every heavy-atom bond. See [Coordinate sets](coordinates.md).
- `--no-phase-angles` uses the raw torsions. By default each torsion that shares its central bond with an earlier torsion is replaced by its offset from that torsion. See [Coordinate sets](coordinates.md).
- `--jacobian` also prints `Mean ln J = ...`, the mean of the log BAT Jacobian over the frames used, on the line before the entropy. It requires `--tree-bonds` and cannot be combined with `--torsions-only`.
- `--legacy` reproduces the original Python implementation bit for bit. See [Legacy mode](legacy.md).
- `--estimator` chooses the estimator: `nn` (the default) is the nearest-neighbor MIE, `qh` is quasi-harmonic, and `bqh` is Boltzmann quasi-harmonic. See [Quasi-harmonic estimators](quasi-harmonic.md). `--mie-order` applies only to `nn`.
- `--histogram-bins N` makes `bqh` estimate each coordinate's first-order entropy from an N-bin histogram instead of nearest neighbors.
- `--periodic-correlation complex|centered` sets how torsions enter the quasi-harmonic correlation term; `complex`, the default, is Harpole and Sharp's `e^{iφ}`. See [Quasi-harmonic estimators](quasi-harmonic.md).
