# Validation

This page records tests of the estimators
against entropies computed without any density estimate, so their accuracy
can be judged on real molecular data rather than only on the closed-form
distributions of the test suite.

## Butane end states

### Reference

Two states of butane differ only in their bonded potential: state 0 has the
torsion terms switched off, state 1 is ordinary butane. Both were simulated
with OpenMM 8.5 in the gas phase (no cutoff) at 320 K with a Langevin
integrator (2 fs, friction 1/ps) and the X-H bonds constrained.

Thermodynamic integration and MBAR over the linear path
U(λ) = (1 - λ) U0 + λ U1, in 21 windows of 10 ns each, give the free energy
difference ΔA. The mean potential energies of the end states give ΔU, and
T ΔS = ΔU - ΔA. Because the masses do not change, the momentum entropy
cancels, and ΔS is the difference in configurational entropy. The metric
term of the constraints is -0.0003 nats, negligible. Three replicates
started from independent frames with independent random seeds:

| replicate | ΔA (kcal/mol) | ΔU (kcal/mol) | S0 - S1 (nats) |
|---|---|---|---|
| 1 | 10.592 ± 0.006 | 1.550 ± 0.051 | 14.219 ± 0.080 |
| 2 | 10.613 ± 0.006 | 1.438 ± 0.051 | 14.429 ± 0.081 |
| 3 | 10.619 ± 0.006 | 1.509 ± 0.059 | 14.325 ± 0.093 |
| mean ± sd | | | 14.32 ± 0.11 |

### Estimates

Each estimate uses the 10,000 frames (one per ps) of the λ = 0 and λ = 1
windows of the same replicate. With the X-H bonds constrained, the default
`CoordinateOptions` give exactly
the 26 unconstrained, nonredundant coordinates of butane: 3 C-C bonds, 12
angles and 11 torsions. The reference is a Cartesian entropy, so every
estimate below includes the change in the mean log BAT Jacobian,
Δ⟨Σ 2 ln b + Σ ln sin θ⟩ over the unconstrained bonds and the angles, about
-0.06 nats; the constrained bonds contribute a constant that cancels.

`S0 - S1` in nats. Mean error and RMSE are against each replicate's own
reference.

| coordinates | estimator | rep 1 | rep 2 | rep 3 | mean error | sd | RMSE |
|---|---|---|---|---|---|---|---|
| | **reference** | 14.219 | 14.429 | 14.325 | | | |
| phase angles | **BQH** | 14.238 | 14.480 | 14.348 | **+0.031** | 0.017 | **0.034** |
| phase angles | BQH, centered torsions | 14.214 | 14.459 | 14.327 | +0.009 | 0.019 | 0.018 |
| phase angles | BQH, 100-bin histogram | 14.165 | 14.340 | 14.226 | -0.081 | 0.024 | 0.083 |
| phase angles | NN MIE-2 | 14.541 | 13.950 | 14.160 | -0.107 | 0.403 | 0.346 |
| phase angles | QH | 13.388 | 13.505 | 13.408 | -0.891 | 0.051 | 0.892 |
| raw torsions | QH | 13.730 | 13.651 | 13.620 | -0.657 | 0.150 | 0.669 |
| raw torsions | NN MIE-2 | 16.315 | 15.446 | 17.879 | +2.222 | 1.273 | 2.453 |
| raw torsions | BQH, 100-bin histogram | 17.132 | 17.331 | 17.167 | +2.886 | 0.038 | 2.886 |
| raw torsions | BQH | 17.189 | 17.446 | 17.147 | +2.936 | 0.102 | 2.938 |
| raw torsions | BQH, centered torsions | 20.175 | 16.012 | 15.787 | +3.000 | 2.560 | 3.657 |

BQH is `estimate_quasi_harmonic_entropy`
with the default options (nearest-neighbor marginals, complex torsions);
"centered torsions" is `PeriodicCorrelation::Centered`
and the histogram is `MarginalEstimator::Histogram`.
NN MIE-2 is `calculate_entropy` at order 2.
"Raw torsions" sets `CoordinateOptions::phase_angles`
to false.

### Split estimator

The [split estimator](split.md) with the topological soft set, on the same trajectories (`S0 - S1` in nats). 96 configurations were scored: soft estimator {MIE-2, joint with either finite-sample constant} × hard estimator {QH, BQH} × hard phase angles {centered, complex} × coupling {per-state, Gaussian e^{iφ}, Gaussian cos/sin, none} × soft set {topology, detected}.

| estimator | rep 1 | rep 2 | rep 3 | mean error | RMSE |
|---|---|---|---|---|---|
| split (MIE-2 soft, QH hard, Gaussian coupling with three harmonics) | 14.296 | 14.565 | 14.482 | +0.123 | 0.128 |
| split, cos/sin coupling | 14.216 | 14.499 | 14.408 | +0.050 | 0.063 |
| split, per-state coupling | 14.246 | 14.529 | 14.440 | +0.081 | 0.089 |
| split, joint soft, BQH hard, no coupling | 14.220 | 14.557 | 14.368 | +0.057 | 0.078 |

On butane the coupling between soft and hard coordinates is small (with cos/sin or `e^{iφ}` soft features, Gaussian and no coupling differ by 0.001 nats), detection chooses the same soft set as topology, and the finite-sample constant cancels in the difference. Centered hard phase angles beat complex ones by about 0.02 nats, and per-state coupling is about 0.03 nats worse than Gaussian coupling, from the bias of a Gaussian fitted within each state. The split estimator improves on an earlier hand-tuned version of the same decomposition (soft torsions by a joint estimate or MIE-2, hard Gaussian within hand-binned states: +0.06 and +0.14 nats), but plain BQH is better on butane, whose single multi-basin torsion BQH's nearest-neighbor marginals already describe. The default's three harmonics, chosen on the molecules below, add about 0.07 nats here, as they do at 300 K.

### What it shows

- **BQH with phase angles agrees with the reference within the reference's
  own uncertainty** (about 0.08 nats per replicate), and follows its
  variation from replicate to replicate. The centered and complex torsion
  treatments are indistinguishable here: with phase angles only four
  torsions remain free, and their correlations are weak.
- **QH misses by 0.9 nats.** It shares BQH's correlation term, so the
  difference is its Gaussian first-order term, which cannot represent the
  multiwell torsions of state 1.
- **NN MIE-2 is right on average but noisy**, ±0.4 nats between replicates:
  the pairwise mutual informations it sums each carry estimator error.
- **Phase angles matter for every estimator.** With the raw torsions, the
  three hydrogen torsions of each methyl group are redundant; every
  estimator is off by 0.7 to 3.7 nats.

These are estimates for one small molecule. They show which estimators work
on it, not an error bar for other systems.

### Reproducing

`examples/butane_validation.rs` computes every estimate in the table from
the end-state trajectories of each replicate:

```bash
cargo run --release --example butane_validation -- \
    rep1/endpoints 14.219 rep2/endpoints 14.429 rep3/endpoints 14.325
```

Each directory holds `0.0.top`, `0.0.nc`, `1.0.top` and `1.0.nc`; the number
after it is that replicate's reference `S0 - S1`.

## Small molecules and peptides at 300 K

Butane has one multi-basin torsion. To set the split estimator's defaults and find where it fails, the same comparison was repeated on 18 small molecules and two peptides in two environments, all at 300 K. The defaults were then checked on 8 more small molecules that played no part in choosing them.

### Systems and end states

- **Small molecules**, GAFF2 with AM1-BCC charges, gas phase: alkanes (butane, hexane, decane, hexadecane), molecules from the FreeSolv set, and others chosen as edge cases: a linear group (1-butyne), an intramolecular hydrogen bond (ethylene glycol), a flat double well (biphenyl), a puckering ring (cyclohexane) and a gas-phase zwitterion (GABA). Held out (‡): 3-amino-1-propanol, 1,2-dimethoxyethane, 2,3-dimethylbutane, ethyl acetate, N-ethylacetamide, 1-pentanol, butylbenzene and alanine dipeptide (Ace-Ala-Nme, with GAFF2 like the others).
- **Peptides**, ff14SB: chignolin (CLN025), a folded β-hairpin, and Ace-(GS)<sub>5</sub>-Nme, a disordered chain, each in vacuum and in GBn2 implicit solvent (mbondi3 radii).
- **State 1** is the molecule with its full force field. **State 0** keeps bonds and angles but has every dihedral term, charge and Lennard-Jones parameter switched off, and in GBn2 the solvation term too. For the molecules marked † state 0 keeps impropers and the torsions about ring bonds, so rings stay planar and only rotatable bonds become free; for the two marked "ring torsions off" the ring is floppy in state 0.

### Reference

Hamiltonian replica exchange over λ windows (neighbor swaps every ps) with MBAR, along a path that scales the torsions linearly and the nonbonded terms with a soft core (Beutler, α = 0.5), on λ<sub>k</sub> = sin²(πk/2(K−1)). Small molecules: 21 windows of 10 ns, 2 fs steps, X-H bonds constrained; butane and hexane used independent windows and thermodynamic integration instead. Peptides: 32 (GS) or 40 (chignolin) windows of 30 ns, 4 fs steps with hydrogen mass repartitioning. T ΔS = ΔU − ΔA as for butane. Three replicates each; the table gives the mean and standard deviation of S0 − S1.

Checks on the references:

- **Replicates** agree to 0.01–0.12 nats for the small molecules and 0.04–0.4 for the GS chain. Chignolin's differ by up to 6 nats: its folding equilibrium is not converged, so its errors below are good to about ±3 nats.
- **Finite differences**: for three molecules simulated identically at 300 and 320 K, −∂ΔA/∂T from the two free energies agreed with (ΔU − ΔA)/T within its error (0.03–0.41 nats against errors of 0.2–0.3), which rules out a problem in the mean energies.
- **Constraints**: four molecules rerun without constraints (1 fs) gave references 0.03–0.07 nats higher. This is within the noise of the comparisons below.

References are therefore good to about ±0.1 nats for the small molecules; differences between estimators smaller than that are not resolved.

### Estimates

Each estimate uses the end-state frames (one per ps) of the same replicate, phase angles and the change in ⟨ln J⟩, as for butane. RMSE in nats over the three replicates against each replicate's own reference; "soft" is the number of topological soft torsions. "Split" is the default (MIE-2 soft, centered QH hard, Gaussian coupling with three harmonics of each soft torsion); "cos/sin" uses only the first harmonic in the coupling (the default before the held-out check), "MIST" replaces MIE-2 for the soft set, and "per-state" replaces the coupling. The best estimate in each row is in bold. "Fails" means no state had enough frames.

| system | soft | reference | QH | BQH | split | split mean error | split, cos/sin | split, MIST | split, per-state |
|---|---|---|---|---|---|---|---|---|---|
| butane | 3 | 14.51 ± 0.03 | 1.03 | 0.11 | 0.09 | +0.08 | **0.04** | 0.09 | 0.06 |
| hexane | 5 | 22.01 ± 0.12 | 1.67 | 0.07 | 0.13 | +0.12 | **0.04** | 0.12 | 0.50 |
| decane | 9 | 12.17 ± 0.11 | 7.43 | 1.95 | 0.73 | -0.34 | 1.30 | **0.66** | fails |
| hexadecane | 15 | 21.22 ± 0.12 | 13.13 | 4.34 | **0.40** | +0.18 | 1.22 | 1.52 | fails |
| 1-butyne | 3 | 0.98 ± 0.02 | 0.69 | 0.11 | 0.05 | +0.01 | 0.07 | **0.04** | 0.07 |
| dichloropropane | 2 | 2.58 ± 0.03 | 1.79 | 0.54 | **0.06** | -0.06 | 0.22 | 0.06 | 0.15 |
| ethylene glycol | 3 | 2.96 ± 0.01 | 2.29 | 1.42 | **0.43** | -0.43 | 0.56 | 0.66 | 0.48 |
| isoflurane | 3 | 4.33 ± 0.08 | 3.21 | 1.50 | **0.14** | -0.13 | 0.29 | 0.21 | 0.27 |
| butylamine | 4 | 4.19 ± 0.05 | 2.77 | 0.55 | **0.06** | -0.05 | 0.30 | 0.08 | 0.18 |
| trimethyl phosphate | 6 | 3.77 ± 0.07 | 2.68 | 0.72 | **0.22** | -0.12 | 0.38 | 0.29 | 0.41 |
| methylpentanol | 6 | 7.48 ± 0.07 | 4.40 | 1.17 | **0.20** | -0.12 | 0.60 | 0.37 | 0.49 |
| trimethoxyethane | 7 | 9.61 ± 0.03 | 8.13 | 3.27 | 0.33 | +0.26 | **0.20** | 1.22 | 0.22 |
| trifluoroethylaniline | 4 | 3.31 ± 0.09 | 1.82 | 1.14 | **0.23** | -0.21 | 0.61 | 0.34 | 0.49 |
| biphenyl † | 1 | 1.40 ± 0.06 | 0.70 | 0.30 | 0.17 | -0.17 | 0.31 | 0.17 | **0.02** |
| N,N-dimethyl-4-nitrobenzamide † | 5 | 6.90 ± 0.11 | 5.80 | 3.50 | 0.84 | -0.84 | 1.63 | 1.03 | **0.34** |
| N,N-dimethyl-4-nitrobenzamide, ring torsions off | 5 | 20.11 ± 0.11 | **0.24** | 2.28 | 3.75 | +3.74 | 2.96 | 3.50 | 4.28 |
| cyclohexane, ring torsions off | 0 | 4.53 ± 0.10 | 0.72 | 6.77 | 0.61 | +0.61 | 0.61 | 0.61 | **0.61** |
| GABA zwitterion | 4 | 9.58 ± 0.06 | 6.58 | 2.73 | **1.00** | -0.99 | 1.81 | 1.78 | 3.67 |
| 3-amino-1-propanol ‡ | 4 | 6.98 ± 0.04 | 4.74 | 2.36 | **0.46** | -0.44 | 0.79 | 1.22 | 0.83 |
| 1,2-dimethoxyethane ‡ | 5 | 4.94 ± 0.04 | 3.78 | 1.21 | **0.41** | -0.41 | 0.66 | 0.59 | 0.52 |
| 2,3-dimethylbutane ‡ | 5 | 6.75 ± 0.06 | 4.78 | 1.22 | **0.43** | -0.40 | 0.92 | 0.47 | 0.86 |
| ethyl acetate ‡ | 4 | 5.16 ± 0.04 | 1.40 | 1.00 | 0.17 | +0.16 | 0.25 | **0.10** | 0.73 |
| N-ethylacetamide ‡ | 4 | 5.96 ± 0.11 | 2.70 | 0.65 | 0.34 | +0.32 | 0.30 | **0.25** | 0.56 |
| 1-pentanol ‡ | 5 | 5.11 ± 0.06 | 3.01 | 0.85 | 0.22 | -0.15 | 0.49 | 0.21 | **0.11** |
| butylbenzene † ‡ | 4 | 4.66 ± 0.04 | 2.98 | 1.07 | **0.22** | -0.22 | 0.57 | 0.28 | 0.28 |
| alanine dipeptide (Ace-Ala-Nme) ‡ | 7 | 14.28 ± 0.02 | 5.19 | 2.06 | 0.75 | +0.60 | 0.55 | **0.08** | 0.94 |
| Ace-(GS)5-Nme, GBn2 | 43 | 104.02 ± 0.04 | 31.09 | 22.35 | **7.04** | -5.55 | 9.38 | 14.16 | fails |
| Ace-(GS)5-Nme, vacuum | 43 | 136.86 ± 0.37 | 50.05 | 40.55 | 57.28 | +56.95 | 53.37 | **28.93** | 63.66 |
| chignolin CLN025, GBn2 | 50 | 194.59 ± 3.00 | 21.45 | 5.30 | 30.34 | +29.46 | 26.98 | **3.61** | fails |
| chignolin CLN025, vacuum | 50 | 231.73 ± 3.31 | 53.44 | 35.00 | 164.76 | +163.98 | 158.41 | **28.13** | 193.72 |

QH and BQH are `estimate_quasi_harmonic_entropy` with complex torsions.

### What it shows

- **On small flexible molecules the split estimator is the best general choice.** On the 23 small molecules with only rotatable bonds free in state 0 (all but the two "ring torsions off" rows and GABA), the default was within 0.5 nats on 20, with an RMSE of 0.38 and 0.84 at worst. It beat BQH on all but hexane. QH misses multi-basin torsions by up to 13 nats.
- **Three harmonics in the coupling held up on new molecules.** On the 15 molecules used to choose it, three harmonics instead of cos/sin cut the RMSE from 0.70 to 0.36 nats. On the 8 held-out ones it cut it from 0.61 to 0.41, improving 6 of the 8. Higher harmonics let a hard coordinate depend on which well a torsion is in, or on cos 3φ as an eclipsed bond angle does. Cos/sin alone left S0 − S1 too low on almost every molecule. It is still better on butane, hexane and trimethoxyethane, by up to 0.13 nats, and on alanine dipeptide, where three harmonics overshoot by 0.6. On each end state alone the coupling is biased downward by about (soft features × hard coordinates) / 2n, as for any Gaussian mutual information. That bias is the same in both end states and cancels in S0 − S1, which was unchanged from 1,250 to 30,000 frames.
- **MIST** (`--soft-estimator mist`) is the soft estimator for folded peptides. It was within 3.6 nats on chignolin in GBn2 (BQH 5.3) and 28 in vacuum (BQH 35). It also roughly halved the error on the GS chain in vacuum, and it nailed alanine dipeptide (0.08). MIE-2 subtracts every pairwise mutual information and so overcounts dependence when many torsions are coupled, which is why the default overestimated the entropy change of chignolin by 30–165 nats. On small molecules MIST was worse than MIE-2 (RMSE 0.61 against 0.38), most of all with hydrogen bonds (aminopropanol 1.2, trimethoxyethane 1.2). On the disordered GS chain in GBn2 it was worse too (14 against 7).
- **Per-state coupling** is better when there are few, well-sampled states (biphenyl 0.02, the ring-kept benzamide 0.34), and worse or impossible with many (hexane, decane and longer).
- **Compact molecules are hardest.** The GABA zwitterion, which folds into a salt bridge, was off by 1.0 nats with only 4 soft torsions.
- **Rings**: with the topological soft set a puckering cyclohexane falls to QH (0.6 nats); BQH was off by 6.8 there. Soft-set detection with MIE-2 was off by 8–10 nats, because the pairwise MIE cannot represent the ring closure; the joint soft estimator was within 0.14. With a benzamide's ring torsions off, plain QH was best.
- **Linear groups** (1-butyne) need nothing special.

### Best estimator for each system

The lowest RMSE in each row of the table above, and the default's where it was not the best.

| system | best | RMSE | default | notes |
|---|---|---|---|---|
| butane, hexane | split, cos/sin | 0.04 | 0.09, 0.13 | one or two multi-basin torsions; BQH also within 0.11 |
| decane | split, MIST | 0.66 | 0.73 | |
| hexadecane | split (default) | 0.40 | | |
| 1-butyne, dichloropropane | split, MIST or default | 0.04–0.06 | 0.05, 0.06 | |
| ethylene glycol, isoflurane, butylamine, trimethyl phosphate, methylpentanol, trifluoroethylaniline | split (default) | 0.06–0.43 | | |
| trimethoxyethane | split, cos/sin | 0.20 | 0.33 | |
| biphenyl †, nitrobenzamide † | split, per-state | 0.02, 0.34 | 0.17, 0.84 | one or two state-defining torsions with well-sampled states |
| nitrobenzamide, ring torsions off | QH | 0.24 | 3.75 | every split variant is 3–4 nats off |
| cyclohexane, ring torsions off | split (any) | 0.61 | | no soft torsions, so the split is QH on centered coordinates; BQH 6.8 |
| GABA zwitterion | split (default) | 1.00 | | next best MIST, 1.8 |
| aminopropanol ‡, dimethoxyethane ‡, 2,3-dimethylbutane ‡, butylbenzene † ‡ | split (default) | 0.22–0.46 | | |
| ethyl acetate ‡, N-ethylacetamide ‡ | split, MIST | 0.10, 0.25 | 0.17, 0.34 | |
| 1-pentanol ‡ | split, per-state | 0.11 | 0.22 | |
| alanine dipeptide ‡ | split, MIST | 0.08 | 0.75 | the default overshoots by 0.6 |
| Ace-(GS)<sub>5</sub>-Nme, GBn2 | split (default) | 7.0 | | MIST 14, BQH 22 |
| Ace-(GS)<sub>5</sub>-Nme, vacuum | split, MIST | 29 | 57 | BQH 41 |
| chignolin, GBn2 | split, MIST | 3.6 | 30 | BQH 5.3 |
| chignolin, vacuum | split, MIST | 28 | 165 | BQH 35 |

The default was the best, or within 0.15 nats of the best, on 22 of the 26 small molecules. The exceptions were biphenyl and the ring-kept nitrobenzamide, where per-state coupling was better; the nitrobenzamide with ring torsions off, where QH was better; and alanine dipeptide, where MIST was better. On the peptides MIST was best in three of four. The exception was the disordered chain in implicit solvent.

These results set the defaults and the warning above 20 soft coordinates. They are estimates for these molecules, not error bars for others.
