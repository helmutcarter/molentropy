# Split estimator

Multi-basin torsions are poorly described by a Gaussian, while bonds, angles and phase angles are close to harmonic, and a Gaussian captures all of their correlations at once. The split estimator (`--estimator split`, `estimate_split_entropy`) uses each where it fits. It divides the coordinates into a **soft** set S and a **hard** set H and writes

S = S<sub>soft</sub> + S<sub>hard</sub> + coupling

- **S<sub>soft</sub>**: the nearest-neighbor MIE over the soft coordinates (`--soft-estimator mie`, the default; `--mie-order` sets the order), or one joint nearest-neighbor estimate over all of them (`joint`, up to 8 soft coordinates).
- **S<sub>hard</sub>**: the quasi-harmonic entropy of the hard coordinates over all frames (`--hard-estimator qh`, the default), or BQH (`bqh`). Hard phase angles are centered on their circular mean (`--periodic-correlation centered`, the split default).
- **Coupling**: the dependence between the two sets, at most 0 up to sampling error. See [Coupling](#coupling).

The output lists each term before the total. Add ⟨ln J⟩ for a Cartesian entropy, as for every estimator.

## Which coordinates are soft

`--soft` (`soft=` in Python) chooses the soft set:

- **`topology`** (the default for files and for a Python `System`): the primary torsion about each rotatable bond, that is, about every bond not on a ring. With phase angles, the default, these are the coordinates of type `torsion`; methyl rotors count. Torsions about ring bonds, and every phase angle, bond and angle, are hard. This needs phase angles, because the raw secondary torsions of a rotor are as multimodal as the primary one.
- **`detect`** (the default for a raw array in Python): every periodic coordinate whose distribution has more than one basin, or is a free rotor (mean resultant length below 0.3), as found by [basin detection](#basin-detection).
- **a list** of coordinate numbers, counted from 1 on the command line (as printed by `--per-coordinate`) and from 0 in Python and Rust.

In Rust, `InternalCoordinates::soft_torsions()` gives the topological set; in Python, `System.soft_torsions`. `symmetric_rotor_torsions` lists the soft torsions about symmetric rotors (one end carries only two or more terminal atoms of one element, as in CH<sub>3</sub>, NH<sub>3</sub><sup>+</sup> and CF<sub>3</sub>), and `ring_torsions` those about ring bonds.

## Coupling

`--coupling` (`coupling=`) chooses how the dependence between the sets is accounted for:

- **`gaussian`** (the default): the Gaussian mutual information between the sets, `½ [ln|C(S∪H)| − ln|C(S)| − ln|C(H)|]`, with soft torsions entering the correlation matrix as cos φ and sin φ (`--soft-features cossin`, the default) or as `e^{iφ}` (`complex`). As two real columns, cos φ and sin φ can correlate with a hard coordinate separately: a bond angle that responds to trans against gauche but not to which gauche correlates with cos φ alone, which the single complex entry partly averages away.
- **`states`**: the hard entropy is computed within each conformational state of the soft torsions, the states' entropies are averaged with their populations, and the coupling is the difference from the hard entropy over all frames. A state is the combination of basins of the state-defining torsions: by default (`--state-torsions auto`) the soft torsions with more than one basin, excluding symmetric rotors, whose basins are equivalent; or a list. States with fewer than `--min-state-frames` frames (300 by default) or than the number of hard coordinates plus 2 are dropped and the others' populations renormalized; the output lists every state and the dropped frames, and warns when more than 5% are dropped. In Python, `state_labels` gives each frame's state directly.
- **`none`**: the sets are treated as independent.

No mixing term is added: the entropy of moving between states is already part of S<sub>soft</sub>.

## Basin detection

`detect_basins` (Rust `basins::detect_basins`, Python `nn_entropy.detect_basins`) finds the basins of a periodic coordinate:

1. Histogram the values over one period (72 bins, fewer for small samples) and smooth with a wrapped Gaussian (standard deviation 2 bins).
2. Take every circular local maximum as a peak.
3. Merge neighboring peaks whose separating minimum is at least half the lower peak.
4. Put the basin boundaries at the remaining minima.
5. Merge basins holding less than 2% of the frames into a neighbor.

The result is deterministic and does not depend on where the period starts.

## Choosing settings

The defaults (topological soft set, MIE-2 soft, centered QH hard, Gaussian coupling with cos/sin features) were chosen on 18 small molecules and 4 peptide systems with MBAR references; see [Validation](validation.md#small-molecules-and-peptides-at-300-k).

- **Small flexible molecules** (up to about 15 soft torsions): on 15 such molecules the defaults were within 0.6 nats on 11 and 1.6 at worst, and beat both QH and BQH on all but biphenyl, where they tied with BQH. (With a benzamide's ring torsions switched off as well, QH was better.) With one multi-basin torsion, as in butane, plain BQH is about as good.
- **Per-state coupling** (`--coupling states`) is more accurate when there are a few well-sampled conformers: about 4 or fewer state-defining torsions and plenty of frames in each state. It cut the error from 0.31 to 0.02 nats on biphenyl and from 1.6 to 0.3 on a nitrobenzamide. With more states it drops frames or fails, and on hexane it was 0.5 nats worse than Gaussian coupling.
- **Folded or compact molecules**: on a folded β-hairpin (chignolin) the split estimator overestimated the entropy loss by 25–160 nats, and plain BQH did better (5 nats in implicit solvent). A zwitterion collapsed into a salt bridge was off by 1.8 nats with only 4 soft torsions. Above 20 soft coordinates the CLI and Python warn (`LARGE_SOFT_SET` in Rust); compare with BQH there.
- **Flexible rings**: torsions about ring bonds are hard under the topological rule, so a puckering ring falls to the Gaussian; on cyclohexane this was within 0.6 nats. Do not combine `--soft detect` with MIE-2 there: detection makes the ring torsions soft, and the pairwise MIE cannot represent the ring closure that couples them (8–10 nats off). The joint soft estimator handled it (0.14 nats).
- **The joint soft estimator** has no truncation error but its bias grows with the number of soft torsions, and it is limited to 8.
- **Detection** can call a rotatable bond stuck in one well hard, and a stiff torsion that samples two wells soft; the topological rule cannot.

`--per-coordinate` and `--mutual-information` do not apply to the split estimator.
