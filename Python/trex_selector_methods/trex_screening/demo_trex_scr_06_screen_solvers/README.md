# Demo 06: Screen-TRex Solver Comparison

## Purpose

Does the underlying T-Rex solver backend change screening performance? Screening thresholds the dummy-based voting proportion `Phi_j`, but `Phi` itself is produced by an inner T-Rex solver. This demo runs two Monte Carlo SNR sweeps that hold the screening machinery fixed while swapping the solver: **TLARS** (T-Rex Least-Angle Regression Solver, the default), **TAFS** (T-Rex Adaptive Forward Selection, `rho_afs = 0.3`), and **TOMP** (T-Rex Orthogonal Matching Pursuit) — first under Ordinary screening only, then crossed with both the Ordinary and Bootstrap-CI thresholds. Mirrors `cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/`.

## Data-generating process

`dgp_iid_snr` (from the shared suite-level `trex_scr_common.py`) — i.i.d. standard-Gaussian design with SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.

- **Both parts**: `n = 300`, `p = 1000` (high-dimensional, `p > n`), random 0-based support of size `s = 10` drawn fresh per trial (matching the C++ demo), coefficients `beta_j = 1`, `SNR in {0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, `NUM_MC = 20` (the C++ demo uses 500; reduced for a quick runtime).

## Methods compared

Three solver backends, configured via `make_scr_trex_ctrl()`:

- **TLARS** — `solver="TLARS"` (default; no extra parameters).
- **TAFS-0.3** — `solver="TAFS"`, `rho_afs=0.3`, with stagnation stop enabled (`tloop_stagnation_stop=True`, `tloop_max_stagnant_steps=5`).
- **TOMP** — `solver="TOMP"`, with the same stagnation stop enabled.

Each solver is paired with a screening rule via `make_screen_ctrl(trex_method="TREX", bootstrap=...)`:

- **Ordinary** — `bootstrap=False` (`Phi_j > 0.5`).
- **Bootstrap-CI** — `bootstrap=True` (library defaults for the bootstrap replicates / CI grid).

Selection is performed by `ts.TRexScreeningSelector`.

## Parts

- **Part 1** — MC SNR sweep, three solvers under Ordinary screening only (3 series): `TLARS (Ordinary)`, `TAFS-0.3 (Ordinary)`, `TOMP (Ordinary)`. Writes `simulation_results/data/d06_screen_trex_solvers_mc_snr_n300_p1000_s10.txt` (summary table) and `.csv` (tidy long format).
- **Part 2** — MC SNR sweep, three solvers x {Ordinary, Bootstrap-CI} (6 series). Writes `simulation_results/data/d06_screen_trex_solver_method_mc_snr_n300_p1000_s10.txt` and `.csv`.

Each series averages FDR / TPR / Estimated FDR over 20 runs per (variant, SNR) point.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_solvers/demo_trex_scr_06_screen_solvers.py
```

The script bootstraps `sys.path` with its own directory and the parent suite directory to import the shared `trex_scr_common.py`; the `RUN_PART_1` / `RUN_PART_2` flags at the top toggle individual parts. Results are statistically comparable to the C++ reference, not bit-identical (different RNG streams, and `NUM_MC = 20` vs. 500 upstream, so the curves are noisier than the C++ figures). The `simulation_results/` folder is empty save for `.gitkeep` until the demo is run.

## Interpretation

- **Solver choice should barely move the Ordinary curves.** Expect TLARS, TAFS-0.3, and TOMP to track each other closely: high FDR with near-zero TPR at very low SNR, converging toward low FDR and TPR near 1.0 as SNR rises. This matches the expectation that the screening decision depends mainly on `Phi_j` rather than on solver-specific path details — differences are within Monte Carlo noise at 20 runs. The Estimated FDR tracks the realized FDR.
- **The 6-series part is where a solver interaction may appear — via the bootstrap step, not the path algorithm.** The Ordinary rows in Part 2 reproduce Part 1. Check whether the Ordinary-vs-Bootstrap gap is uniform across solvers, or whether one solver (e.g. TOMP) interacts differently with the bootstrap thresholding step than the LARS-style backends do.
- `TAFS-0.3`'s extra regularization (`rho_afs = 0.3`) is the most likely place a solver-specific difference could surface, since it alters the T-loop's adaptive forward-stagewise behavior rather than just the path algorithm; at 20 runs any effect is likely small relative to noise.

**Last updated**: 2026-07-08
