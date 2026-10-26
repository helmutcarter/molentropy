# Interpretation and limitations

- Results are differential internal-coordinate entropies in natural-log units (equivalently, units of `k_B`), not kcal mol⁻¹ K⁻¹.
- The reported value does not include a BAT-to-Cartesian Jacobian correction, momentum entropy, or a standard-state term, and must not be described as an absolute thermodynamic entropy. `--jacobian` reports the Jacobian term separately, for the tree bond set, but does not add it to the total.
- MD frames are generally correlated, which biases nearest-neighbor estimates; see [Correlated frames](correlated-frames.md).
- The quasi-harmonic correlation term is a Gaussian approximation: it captures correlations in the second moments of the coordinates (for torsions, of `e^{iφ}`), not nonlinear dependence. See [Quasi-harmonic estimators](quasi-harmonic.md).
- The MIE is truncated; see [Nearest-neighbor MIE](nn-mie.md#notes).
- The split estimator's topological soft set misses rotors that the bond graph cannot identify as symmetric, such as t-butyl (whose methyls are not terminal atoms), so their basins can define states; detection can call a rotatable bond stuck in one well hard. Its per-state coupling fragments the frames and fits a Gaussian within each state, which is biased for small states. See [Split estimator](split.md).
- The split estimator was validated on small, flexible molecules, where it was within 0.6 nats on 11 of 15 and 1.6 nats at worst, with up to 15 soft torsions. On a folded peptide it overestimated the entropy change by tens of nats, where BQH did far better, and compact molecules such as a gas-phase zwitterion are less accurate even when small. It warns above 20 soft coordinates. See [Validation](validation.md#small-molecules-and-peptides-at-300-k).
