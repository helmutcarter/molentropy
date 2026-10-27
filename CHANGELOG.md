# Changelog

All notable changes to molentropy are recorded here. Versions follow
[semantic versioning](https://semver.org/); before 1.0, a minor version may
change the API or the default estimates.

## [0.1.0] - 2026-10-26

First release, under the name molentropy (previously nn-entropy).

- Reads Amber topologies (`.parm7`) and NetCDF trajectories (`.nc`) and converts
  every frame to bond-angle-torsion coordinates, with phase angles and an
  optional Jacobian term.
- Estimators:
  - nearest-neighbor mutual information expansion (NN MIE) to fourth order, and
    pairwise mutual information;
  - mutual information spanning tree (MIST);
  - joint nearest-neighbor entropy up to eight dimensions;
  - quasi-harmonic (QH) and Boltzmann quasi-harmonic (BQH);
  - a split estimator: soft torsions by MIE, MIST or the joint estimate, stiff
    coordinates by QH or BQH, and a Gaussian (Fourier-harmonic) or per-state
    coupling term.
- Per-coordinate entropies for every estimator except the joint estimate.
- Theiler windows for time-correlated frames, and exact finite-sample constants.
- Results are bit-for-bit reproducible across runs and thread counts.
- A command-line program, a Rust library and Python bindings.
- Validation against HREX/MBAR reference entropies of small molecules and
  peptides at 300 K.

[0.1.0]: https://github.com/helmutcarter/molentropy/releases/tag/v0.1.0
