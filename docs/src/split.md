# Split estimator

Multi-basin torsions are poorly described by a Gaussian, while bonds, angles and phase angles are close to harmonic, and a Gaussian captures all of their correlations at once. The split estimator (`--estimator split`, `estimate_split_entropy`) uses each where it fits. It divides the coordinates into a **soft** set S and a **hard** set H and writes

S = S<sub>soft</sub> + S<sub>hard</sub> + coupling

- **S<sub>soft</sub>**: the nearest-neighbor MIE over the soft coordinates (`--soft-estimator mie`, the default; `--mie-order` sets the order), the mutual information spanning tree (`mist`), or one joint nearest-neighbor estimate over all of them (`joint`, up to 8 soft coordinates). MIST uses the same single and pairwise nearest-neighbor entropies as MIE-2, but subtracts only the pairwise mutual information along the maximum spanning tree of the soft coordinates: one dependence per coordinate, the strongest that keeps the tree connected. With exact terms it is an upper bound on the joint entropy, while MIE-2 can fall far below it when many torsions are coupled.
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

- **`gaussian`** (the default): the Gaussian mutual information between the sets, `½ [ln|C(S∪H)| − ln|C(S)| − ln|C(H)|]`. Soft torsions enter the correlation matrix as cos kφ and sin kφ for k = 1 to 3 (`--soft-features harmonics`, the default; `--harmonics N` sets the number), as cos φ and sin φ only (`cossin`, the same as `--harmonics 1`), or as `e^{iφ}` (`complex`). As separate real columns, cos φ and sin φ can correlate with a hard coordinate separately: a bond angle that responds to trans against gauche but not to which gauche correlates with cos φ alone, which the single complex entry partly averages away. Higher harmonics catch dependences that are not linear in cos φ and sin φ, such as an angle that follows cos 3φ, strained at every eclipsed position. Each harmonic adds two columns per soft torsion, and the frames must outnumber all the columns.
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

The defaults (topological soft set, MIE-2 soft, centered QH hard, Gaussian coupling with three harmonics) were chosen on 18 small molecules and 4 peptide systems with MBAR references, and checked on 8 more small molecules; see [Validation](validation.md#small-molecules-and-peptides-at-300-k), which also lists the best estimator for each system.

### Why these defaults

Each default is the setting with the lowest error on small flexible molecules that also works on every system. A simpler setting was preferred when errors were within about 0.05 nats. The numbers are RMSEs in nats over the 15 molecules used to choose the defaults, and over the 8 held out.

- **Topological soft set**, not detection. Detection was a little better on the first 15 molecules (0.59 against 0.70 with cos/sin coupling) but worse on the held-out 8 (0.66 against 0.61). It also makes ring torsions soft, which put a puckering cyclohexane 8–10 nats off with MIE-2. Topology is deterministic, needs no sampling to classify a torsion, and treats a rotatable bond the same whether or not it crossed a barrier in the trajectory.
- **MIE-2 soft entropy**, not MIST or the joint estimate. MIE-2 had an RMSE of 0.38 on the 23 small molecules, against 0.61 for MIST. The joint estimate was no better (0.69 on the 13 molecules where it ran, against 0.70), its bias grows with the number of soft torsions, and it stops at 8, so it failed on decane and longer. MIST is an option because it is far better on folded peptides, where MIE-2 overcounts the dependence between many coupled torsions; see below.
- **QH for the hard set**, with phase angles centered on their circular mean. BQH hard entropies were within 0.02 nats of QH on the first 15 molecules and worse on the held-out 8 (0.71 against 0.61). QH is also faster, and the hard coordinates are close to Gaussian by construction. On butane, centered phase angles beat complex ones by 0.02 nats; on the 15 molecules the difference was 0.01.
- **Gaussian coupling**, not per-state. No coupling at all left S0 − S1 too low (1.16 and 1.01). Per-state coupling was the most accurate on the first 15 molecules where it ran (0.33), but it fails when no state has enough frames (decane, hexadecane, both peptides in GBn2), and it was worse on the held-out 8 (0.66). Gaussian coupling always runs.
- **Three harmonics** of each soft torsion in the coupling. The error fell with each harmonic: 0.70, 0.55 and 0.36 for one, two and three on the first 15 molecules, and 0.61, 0.54 and 0.41 on the held-out 8. One-hot basin indicators added to cos/sin did about as well as two harmonics (0.55 and 0.56). More harmonics were not tested; each adds two correlation columns per soft torsion and needs more frames. Three harmonics are slightly worse than cos/sin on butane and hexane (by 0.05–0.09 nats), where cos/sin was already within 0.04.

- **Small flexible molecules** (up to about 15 soft torsions): on 23 such molecules the defaults were within 0.5 nats on 20 and 0.84 at worst, and beat both QH and BQH on all but hexane, where BQH was 0.06 nats closer. (With a benzamide's ring torsions switched off as well, QH was better.) With one multi-basin torsion, as in butane, plain BQH is about as good.
- **Per-state coupling** (`--coupling states`) is more accurate when there are a few well-sampled conformers: about 4 or fewer state-defining torsions and plenty of frames in each state. It cut the error from 0.31 to 0.02 nats on biphenyl and from 1.6 to 0.3 on a nitrobenzamide. With more states it drops frames or fails, and on hexane it was 0.5 nats worse than Gaussian coupling.
- **Folded peptides**: use `--soft-estimator mist`. With the MIE-2 default, the split estimator overestimated the entropy loss of a folded β-hairpin (chignolin) by 30–165 nats, because MIE-2 subtracts every pairwise dependence among its 50 soft torsions. With MIST it was within 3.6 nats in implicit solvent and 28 in vacuum, closer than BQH (5.3 and 35). On small molecules and on a disordered chain in implicit solvent, MIST was worse than the default. Above 20 soft coordinates the CLI and Python warn (`LARGE_SOFT_SET` in Rust); compare MIST with BQH there.
- **Compact small molecules**: a zwitterion collapsed into a salt bridge was off by 1.0 nats with only 4 soft torsions.
- **Flexible rings**: torsions about ring bonds are hard under the topological rule, so a puckering ring falls to the Gaussian; on cyclohexane this was within 0.6 nats. Do not combine `--soft detect` with MIE-2 there: detection makes the ring torsions soft, and the pairwise MIE cannot represent the ring closure that couples them (8–10 nats off). The joint soft estimator handled it (0.14 nats).
- **The joint soft estimator** has no truncation error but its bias grows with the number of soft torsions, and it is limited to 8.
- **Detection** can call a rotatable bond stuck in one well hard, and a stiff torsion that samples two wells soft; the topological rule cannot.

`--per-coordinate` and `--mutual-information` do not apply to the split estimator.
