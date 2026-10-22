# Interpretation and limitations

- Results are differential internal-coordinate entropies in natural-log units (equivalently, units of `k_B`), not kcal mol⁻¹ K⁻¹.
- The reported value does not include a BAT-to-Cartesian Jacobian correction, momentum entropy, or a standard-state term, and must not be described as an absolute thermodynamic entropy. `--jacobian` reports the Jacobian term separately, for the tree bond set, but does not add it to the total.
- MD frames are generally correlated, which biases nearest-neighbor estimates; see [Correlated frames](correlated-frames.md).
- File-based results differ from the historical Python implementation because torsions are now periodic and phase angles replace redundant torsions. Its asymptotic finite-sample convention is kept as the default. `--legacy` reproduces the historical results exactly.
- The quasi-harmonic correlation term is a Gaussian approximation: it captures correlations in the second moments of the coordinates (for torsions, of `e^{iφ}`), not nonlinear dependence. See [Quasi-harmonic estimators](quasi-harmonic.md).
- The MIE is truncated; see [Nearest-neighbor MIE](nn-mie.md#notes).
