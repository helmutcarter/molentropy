# Finite-sampling bias

Every estimator in the crate is biased at a finite number of frames:

- **Nearest-neighbor terms** (NN MIE, MIST, the joint estimate, BQH's marginals): the bias grows with the dimension of each term, so it is small for MIE-2's one- and two-dimensional terms and larger for the joint estimate over many coordinates.
- **Gaussian terms** (QH, the correlation term of BQH, the split estimator's hard set and coupling): a mutual information between p and q correlated columns estimated from n frames is biased downward by about pq/2n. For the split estimator's coupling term on chignolin (50 soft torsions with three harmonics each, and about 375 hard coordinates), this was 1.9 nats at 30,000 frames and 15 nats at 3,750.

The entropy of one ensemble is therefore rarely converged in the number of frames. On chignolin with MIST soft entropies, each end state's estimate rose by about 57 nats between 1,875 and 30,000 frames.

## It cancels in differences

What is usually reported is a difference between two ensembles of the same molecule, and there most of the bias cancels. It depends mainly on the number of frames and coordinates, which are the same in both, and much less on the ensemble. In the chignolin example, the two end states' 57 nats of drift left their difference changed by 1.4 nats. In the coupling term alone, the difference between end states stayed within 0.04 nats from 3,750 to 30,000 frames for chignolin in GBn2 (0.65 in vacuum), and within 0.04 nats from 1,250 to 10,000 frames for the small molecules tested.

The table gives the error in S0 − S1 against the HREX/MBAR references of the [validation](validation.md#small-molecules-and-peptides-at-300-k), from a sixteenth of the frames and from all of them (means of three replicates, nats). Most entries barely move: what remains with all frames is the estimator's own error, not too few frames.

| system | frames | split (default) | split, MIST | BQH | QH |
|---|---|---|---|---|---|
| butane | 625 → 10,000 | +0.16 → +0.08 | +0.14 → +0.08 | −0.12 → +0.00 | −1.09 → −1.03 |
| hexane | 625 → 10,000 | +0.83 → +0.12 | +0.31 → +0.12 | +0.01 → +0.06 | −1.69 → −1.67 |
| 1-pentanol | 625 → 10,000 | +0.54 → −0.15 | +0.01 → −0.20 | −0.82 → −0.84 | −3.30 → −3.01 |
| methylpentanol | 625 → 10,000 | −0.74 → −0.13 | −0.54 → −0.36 | −1.22 → −1.16 | −4.52 → −4.40 |
| 3-amino-1-propanol | 625 → 10,000 | −0.21 → −0.44 | −1.20 → −1.22 | −2.46 → −2.36 | −4.78 → −4.74 |
| biphenyl | 625 → 10,000 | −0.08 → −0.17 | −0.08 → −0.17 | −0.02 → −0.29 | −0.50 → −0.70 |
| alanine dipeptide | 625 → 10,000 | +0.03 → +0.60 | +0.28 → +0.04 | −1.43 → −2.05 | −4.89 → −5.19 |
| hexadecane | 625 → 10,000 | **+4.13 → +0.18** | −1.77 → −1.49 | −4.46 → −4.34 | −13.56 → −13.13 |
| Ace-(GS)<sub>5</sub>-Nme, GBn2 | 1,875 → 30,000 | −6.37 → −5.55 | −15.36 → −14.15 | −22.59 → −22.35 | −31.13 → −31.07 |
| Ace-(GS)<sub>5</sub>-Nme, vacuum | 1,875 → 30,000 | +53.83 → +56.95 | −30.76 → −28.84 | −40.83 → −40.53 | −50.36 → −50.02 |
| chignolin, GBn2 | 1,875 → 30,000 | **+44.52 → +29.46** | −4.77 → −3.32 | −4.56 → −5.04 | −21.39 → −21.41 |
| chignolin, vacuum | 1,875 → 30,000 | +163.39 → +163.98 | −29.99 → −28.08 | −34.66 → −34.90 | −53.75 → −53.38 |

## Where it does not cancel

Cancellation needs the bias to be about the same in both ensembles. It was not for the default split estimator on molecules with many coupled soft torsions, whose MIE-2 soft entropy subtracts a pairwise mutual information for every pair of torsions:

- **Hexadecane** (15 soft torsions): the error fell from +4.13 to +0.18 nats between 625 and 10,000 frames. The two end states drifted by 10.4 and 14.4 nats, so 3.9 nats did not cancel.
- **Chignolin in GBn2** (50 soft torsions): the end states drifted by 50 and 65 nats, so 15 nats did not cancel, and the error jumped between strides (+44.5, +27.4, +32.8, +28.4, +29.5 from 1,875 to 30,000 frames) rather than settling.
- **Alanine dipeptide**: the default's overshoot grew from +0.03 to +0.60 nats with more frames. With fewer frames it was right only by chance.

With MIST soft entropies the same molecules converged: the difference moved by 0.3 nats on hexadecane and 1.4 on chignolin. Cancellation also needs the two ensembles to have the same coordinates and about the same number of frames. A Gaussian term's bias scales as 1/n, so comparing 5,000 frames with 50,000 leaves most of the smaller run's bias in the difference.

## Checking a result

To check a difference, estimate both ensembles again on every 2nd and every 4th frame (`--stride 2`, `--stride 4`, or `stride=` in `load_system`) and compare the differences, not the separate estimates. Thinning keeps the span of each trajectory, so the conformations sampled stay the same and only the number of frames changes. A difference that moves steadily toward a value is converging; if it moves by less than the differences you want to resolve, it has converged for that purpose. If it changes direction from one stride to the next, as the default split estimator did on chignolin, the estimate is dominated by noise or by an estimator that does not suit the system.

On densely saved trajectories, thinning also moves nearest neighbors apart in time, which changes the bias from [correlated frames](correlated-frames.md). Use a Theiler window, scaled to the same span in time (`W / stride` frames), or check that the frames are decorrelated, before reading a change between strides as finite-sample bias.

None of this detects states the trajectory never visited. A peptide that stays folded for the whole run gives the same answer at every stride. Only independent replicates or enhanced sampling reveal that.
