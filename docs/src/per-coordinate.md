# Per-coordinate entropies

`--per-coordinate` on the command line, and the `estimate_coordinate_*` functions and `per_coordinate=True` in Python and Rust, divide an entropy among the coordinates. A joint entropy belongs to the whole set, so any such division is a convention. The conventions below share three properties:

- the shares sum to the total, up to rounding;
- they do not depend on the order of the coordinates;
- a coordinate independent of all the others keeps its own entropy.

Use them to compare the same coordinate between two ensembles under the same estimator: which bonds, angles or torsions carry an entropy change. Shares from different estimators divide different quantities and are not comparable.

## By estimator

| estimator | each coordinate's share | where |
|---|---|---|
| NN MIE, order 1 | its marginal entropy | `estimate_coordinate_entropy` |
| NN MIE, orders 2 to 4 | its marginal entropy, plus an equal part of every interaction term it is part of | `estimate_coordinate_entropy` |
| MIST | its marginal entropy, minus half the mutual information of each spanning-tree link it is part of | `estimate_coordinate_mist_entropy` (Rust); the split estimator's soft set |
| QH, BQH | its first-order term, plus ½ (log C)<sub>ii</sub> | `estimate_coordinate_quasi_harmonic_entropy`; see [Quasi-harmonic estimators](quasi-harmonic.md#per-coordinate-entropies) |
| split | its share of the soft or hard entropy as above, plus its share of the coupling | `estimate_coordinate_split_entropy` (Rust), `estimate_split_entropy(..., per_coordinate=True)` (Python) |
| joint NN | not available | |

### The MIE at higher orders

The MIE truncated at order k is a sum of interaction terms: the marginal entropies, minus the pairwise mutual informations, plus the three-way co-informations, and so on up to k coordinates at a time. Each term is split evenly among the coordinates it involves. At order 2 this splits each mutual information evenly between its two coordinates. At every order it is the Shapley value of the truncated expansion, the unique split that treats all coordinates alike and gives an independent coordinate nothing from the others. The shares cost no more than the total. They are computed from the same joint entropies, summed per coordinate.

### The split estimator

- **Soft coordinates** take their shares of the soft entropy: MIE at its order, or MIST. The joint soft estimate has no split, so `--per-coordinate` refuses it.
- **Hard coordinates** take their shares of the QH or BQH entropy.
- **Gaussian coupling** `½ [ln|C(S∪H)| − ln|C(S)| − ln|C(H)|]` is shared like the QH correlation term. A soft coordinate receives `½ [(log C(S∪H))ᵢᵢ − (log C(S))ᵢᵢ]`, a hard one the same with C(H). A soft torsion's share is summed over its columns, two per harmonic. Two coordinates coupled only to each other share their coupling evenly. Unlike a QH correlation share, one coordinate's coupling share can be positive; the sum cannot.
- **Per-state coupling** gives each hard coordinate its population-weighted share of the hard entropy within the states, minus its share over all frames. Soft coordinates get no coupling share, because conditioning on the states changes only the hard entropy.

The command line prints the split shares as `coordinate`, `type`, `atoms`, `entropy`, `set` (`soft` or `hard`), `own` and `coupling`, where `entropy` = `own` + `coupling`.

The shares come from an eigendecomposition, which is stricter about nearly singular matrices than the Cholesky factor behind the totals. Soft torsions that barely move make their harmonics nearly collinear (sin kφ is about kφ for small φ), so a coupling that the total can still compute may be refused for shares. Use the cos/sin features (`--soft-features cossin`) or fewer harmonics there. The peptides and small molecules of the [validation](validation.md) had no such problem.
