# PCA estimator (experimental)

> **Experimental.** This estimator is in molentropy so it can be compared with the others on the same benchmarks, and it was far less accurate than them (below). Do not use it for results.

`--estimator pca` (Python `estimate_pca_entropy`, Rust `pca::estimate_pca_entropy`) follows Numata, Wan and Knapp (*Genome Informatics* 18, 192, 2007). Unlike every other estimator here, it works in mass-weighted Cartesian coordinates rather than bond-angle-torsion coordinates:

1. **Superposition.** Every frame of the solute is superposed on the average structure (mass-weighted, iterated until the average stops moving), which removes overall translation and rotation.
2. **Principal components.** The covariance of the mass-weighted coordinates `√m x` is diagonalized. Its `3N − 6` largest eigenvalues `λ_i` (amu Å²) are the variances of the principal components; the six smallest belong to the removed rigid-body motion.
3. **Quasi-harmonic entropy.** Each component is treated as a quantum harmonic oscillator of frequency `ω_i = √(k_B T / λ_i)`. With `α_i = ħω_i / k_B T`, its entropy is `α_i / (e^{α_i} − 1) − ln(1 − e^{−α_i})`. The sum over components is the quasi-harmonic upper bound of Andricioaei and Karplus.
4. **Corrections for the classical modes** (`α_i < 1`). The trajectory is projected on them, and two corrections are subtracted:
   - *Anharmonicity:* each mode's Gaussian entropy `½ ln(2πe λ_i)` minus its k-nearest-neighbor entropy. A correction below 0.007 of the mode's quantum entropy counts as 0.
   - *Pairwise dependence:* the mutual information of every pair of classical modes, from one- and two-dimensional k-nearest-neighbor entropies, with negative values counted as 0.

The k-nearest-neighbor entropy uses the k-th neighbor (`--neighbors`, default 4, as in the paper): `ψ(n) − ψ(k) + ln V_d + (d/n) Σ ln R_k`. The temperature (`--temperature`, default 300 K) must be the simulation's.

The result is an absolute entropy of the vibrations in units of k_B, without the momentum part. Only differences between ensembles of the same molecule compare with the other estimators.

## How it compared

These are the same end states as on the [Validation](validation.md) page: errors in S0 − S1 against the HREX/MBAR or TI references at 300 K, as the mean over three replicates, in nats.

| systems | PCA | PCA, quasi-harmonic term only | split (default) | split, MIST |
|---|---|---|---|---|
| 24 small molecules (RMSE) | 7.19 | 6.13 | 0.85 | 0.99 |
| Ace-(GS)5-Nme, GBn2 | +6.4 | +19.1 | −5.6 | −14.2 |
| Ace-(GS)5-Nme, vacuum | +306 | +74 | +57 | −29 |
| chignolin, GBn2 | +175 | +140 | +29 | −3.3 |
| chignolin, vacuum | +278 | +181 | +164 | −28 |
| Trp-cage, GBn2 | +58 | +378 | −3.3 | −16.7 |

- **The small molecules were 6–8 times worse than with the split estimator.**
  - Butane and hexane were 6–7 nats low and alanine dipeptide 8 nats high.
  - A torsion that jumps between wells moves atoms along arcs. In Cartesian coordinates that is neither Gaussian nor a single principal component, and the corrections for single modes and pairs recover only part of it.
- **The corrections often made things worse.** On 16 of the 24 small molecules, the error was larger after the anharmonic and pairwise corrections than with the quasi-harmonic term alone. Like MIE-2 on torsions, the sum of pairwise mutual information counts the same dependence several times. On the GS chain in vacuum, the corrections took 232 nats more from state 1 than from state 0, moving the error from +74 to +306.
- **On the peptides it was unreliable in both directions.** Its one good result, the GS chain in GBn2 (+6.4), came from a large overestimate (+19) and a large correction (−13) nearly cancelling.
- **The method was published with no thermodynamic reference.** It was compared only with normal modes of a minimized structure. These benchmarks are its first test against reference entropies here.

Other caveats:
- Superposition is ill-defined for a floppy chain, such as the ideal chain of state 0.
- The rigid-body rotation it removes depends on the structure, which a difference between ensembles of different shapes does not cancel.
- The quantum treatment of stiff modes differs from the classical references, though that difference mostly cancels between states.

Settings: the 24 small molecules are those of the validation table that have three replicates of 10,000 frames. The peptides use 30,000 frames.
