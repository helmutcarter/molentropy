# Correlated frames and the Theiler window

The Kozachenko-Leonenko estimator assumes independent samples. When frames are saved faster than a coordinate decorrelates, a frame's nearest neighbor is usually the frame before or after it. That distance measures how far the system moves in one output interval, not the density of states, so entropies come out low. The error grows with the number of dimensions, so it does not cancel in the MIE: oversampled independent coordinates show spurious positive mutual information.

A Theiler window of `W` frames keeps every frame as a query point but takes its nearest neighbor only from frames more than `W` away in time. Each frame's finite-sample term uses the number of frames it can actually compete against, so the window does not bias the estimate for independent samples. Unlike `--stride`, it discards no data. Frames are counted after `--start` and `--stride`, and the series must have at least `2W + 2` frames.

What the tests show (`tests/theiler.rs`):
- On an AR(1) series with lag-one correlation 0.99, the mutual information of two independent coordinates is +0.044 at `W = 0` and about 0 at `W = 5`.
- On mildly correlated series (lag-one correlation 0.95), any `W` of a few frames or more removes the bias and leaves the estimate on a plateau. RMSE is lower than with any stride.
- On strongly oversampled series the window over-corrects as `W` grows. The frames that remain eligible come in clusters, one per visit to a region, which lengthens nearest-neighbor distances. Mutual information then turns negative, about -0.04 at `W >= 50` in the example above. A window is not a substitute for adequate sampling.

Outside the test suite, on 100 ns production runs of three small solvated ligands saved every 2 ps (50,000 frames), the order-2 total changes by at most about 0.4 k_B across `W = 0` to `250`, while striding moves it by 10-20 k_B. At that spacing the frames are already effectively decorrelated, so the plain estimator is fine there.

Use a window when output is dense compared with the coordinates' correlation times. Choose `W` from a sweep: plot the entropy against `W` and take a value where the curve has stopped moving for the fast coordinates. Prefer a window to a stride when the goal is only to remove temporal neighbors.
