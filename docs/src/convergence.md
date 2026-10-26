# Convergence check

Every estimator is biased at a finite number of frames. The nearest-neighbor terms are biased by an amount that grows with their dimension. A Gaussian mutual information between p and q correlated columns, such as the split estimator's coupling term, is biased by about −pq/2n. For chignolin's 50 soft torsions this was 1.9 nats at 30,000 frames. `--convergence` (`check_convergence` in Python and Rust) shows how far an estimate still is from its large-sample limit.

## What it does

It repeats the estimate on every 4th, every 2nd and every frame. Thinning keeps the span of the trajectory, and so the conformations it sampled, and changes only how many frames describe them. A change between the three estimates is therefore finite-sample bias of the estimator, not drift in the sampling. It costs about 1.75 times the estimate alone.

```text
stride	frames	entropy
4	2500	...
2	5000	...
1	10000	...
Change from 5000 to 10000 frames = ...
Extrapolated entropy = ...
Total entropy = ...
```

- **Change**: the estimate from all frames minus that from every 2nd frame. If the bias falls as 1/n, the estimate from all frames is still off by about this much in the same direction. The estimate has converged for your purpose when the change is small next to the differences you want to resolve.
- **Extrapolated entropy**: the intercept of a least-squares line through the three estimates against 1/n. It is the limit for infinitely many frames if the bias falls as 1/n. That holds for Gaussian terms, and roughly for the one- and two-dimensional nearest-neighbor terms of MIE-2 and MIST. A bias that falls more slowly, such as that of the joint nearest-neighbor estimate over many coordinates, is only partly removed. On independent soft and hard sets, where the true coupling is 0, the extrapolation removed the coupling's bias on average. Sampling noise left a residual of up to 28% of it (`tests/convergence.rs`).

If the estimate does not move the same way at both steps, rising from n/4 to n/2 frames and falling from n/2 to n or the reverse, the CLI warns (`steady` is false in Python; `Convergence::steady` in Rust). Either sampling noise is larger than the bias or the bias does not fall steadily. The change is then a rough size of the remaining error, and the extrapolation means little. On chignolin's folded state, the default split estimator rose 8 nats from 7,500 to 15,000 frames, then fell 2 nats to 30,000. With MIST soft entropies it rose steadily (6.7 then 3.3 nats), and S0 − S1 moved by less than 1.3 nats across the three strides.

The split estimator keeps the soft set of the full run at every stride, so detection cannot change the coordinates between the three estimates. A Theiler window of W frames becomes ceil(W/stride) frames on the thinned data, so it spans the same time.

## Differences between ensembles

For a difference between two ensembles with the same coordinates and the same number of frames, much of the bias cancels: in the 300 K validation, the coupling term's bias was the same in both end states, and S0 − S1 did not change from 1,250 to 30,000 frames. Run the check on both ensembles and compare the changes, or difference the extrapolated values.

## What it cannot detect

- **States the trajectory never visited.** If a peptide stays folded for the whole run, every thinned estimate agrees and the answer is still wrong. Only independent replicates or enhanced sampling reveal that.
- **Drift in the sampling.** Thinning keeps the whole span, so it does not compare early and late parts of the trajectory. To check for drift, estimate the first and second halves separately (`--stop`, `--start`).
- **Correlated frames.** When frames are saved faster than the coordinates decorrelate, thinning also moves nearest neighbors apart in time. The change then mixes finite-sample bias with the bias of temporal neighbors. Use a Theiler window, or check that the frames are decorrelated, first; see [Correlated frames](correlated-frames.md).
- Every 4th frame must still give a usable estimate: the split and quasi-harmonic estimators need more frames than correlation-matrix columns, and per-state coupling needs enough frames per state.

## Python and Rust

In Python, `check_convergence(estimate, data)` calls any estimator on the thinned frames:

```python
system = nn_entropy.load_system(top, traj)
data = np.array(system.data)
check = nn_entropy.check_convergence(
    lambda d: nn_entropy.estimate_split_entropy(
        d, soft=system.soft_torsions, periods=system.periods
    ),
    data,
)
check["frames"], check["entropies"], check["change"], check["extrapolated"], check["steady"]
```

`estimate` takes an array of coordinates by frames and returns a number or a dict with a `"total"`. `data` is such an array or a `System`. A Theiler window counts frames of the array `estimate` receives, whose number of frames is `d.shape[1]`.

In Rust, `convergence::check_convergence(data, frames_end, estimate)` calls `estimate(thinned, frames, stride)` and returns a `Convergence` with the three `points`, the `change` and the `extrapolated` value; `steady()` says whether the estimate moved the same way at both steps. `convergence::extrapolate` and `convergence::subsample` are the pieces it uses.
