# Demo 06: Screen-TRex Solver Comparison

## Purpose

Does the underlying T-Rex solver backend change screening performance? Screening thresholds the dummy-based voting proportion `Phi_j`, but `Phi` itself is produced by an inner T-Rex solver. This demo runs two Monte Carlo SNR sweeps that hold the screening machinery fixed while swapping the solver: **TLARS** (T-Rex Least-Angle Regression Solver, the default), **TAFS** (T-Rex Adaptive Forward Selection, `rho_afs = 0.3`), and **TOMP** (T-Rex Orthogonal Matching Pursuit) — first under Ordinary screening only, then crossed with both the Ordinary and Bootstrap-CI thresholds. Mirrors `cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/`.

## Data-generating process

`dgp_iid_snr` — i.i.d. standard-Gaussian design with SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.

- **Both parts**: `n = 300`, `p = 1000` (high-dimensional, `p > n`), random 1-based support of size `s = 10` drawn fresh per trial (the C++ demo uses 0-based supports), coefficients `beta_j = 1`, `SNR in {0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, `num_MC = 20` (the C++ demo uses 500; reduced for a ~1–2 minute runtime).

## Methods compared

Three solver backends, configured via `make_scr_trex_ctrl()`:

- **TLARS** — `solver = "TLARS"` (default; no extra parameters).
- **TAFS-0.3** — `solver = "TAFS"`, `rho_afs = 0.3`, with stagnation stop enabled (`tloop_stagnation_stop = TRUE`, `tloop_max_stagnant_steps = 5`).
- **TOMP** — `solver = "TOMP"`, with the same stagnation stop enabled.

Each solver is paired with a screening rule via `trex_screen_control(trex_method = "TREX", use_bootstrap_CI = ...)`:

- **Ordinary** — `use_bootstrap_CI = FALSE` (`Phi_j > 0.5`).
- **Bootstrap-CI** — `use_bootstrap_CI = TRUE` (library defaults for the bootstrap replicates / CI grid).

Selection is performed by `TRexScreeningSelector`.

## Parts

- **Part 1** — MC SNR sweep, three solvers under Ordinary screening only (3 series): `TLARS (Ordinary)`, `TAFS-0.3 (Ordinary)`, `TOMP (Ordinary)`. Writes `simulation_results/data/d06_screen_trex_solvers_mc_snr_n300_p1000_s10.txt` (summary table) and `.csv` (tidy long format).
- **Part 2** — MC SNR sweep, three solvers x {Ordinary, Bootstrap-CI} (6 series). Writes `simulation_results/data/d06_screen_trex_solver_method_mc_snr_n300_p1000_s10.txt` and `.csv`.

Each series averages FDR / TPR / Estimated FDR over 20 runs per (variant, SNR) point.

## Running

```bash
Rscript R/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_solvers/demo_trex_scr_06_screen_solvers.R
```

Runs from any working directory; the `run_part_1` / `run_part_2` flags at the top toggle individual parts; the script sources the shared `../trex_scr_sim_utils.R`. Results are statistically comparable to the C++ reference, not bit-identical (different RNG streams, and `num_MC = 20` vs. 500 upstream, so the curves are noisier than the C++ figures).

## Interpretation

- **Solver choice barely moves the Ordinary curves.** In the recorded Part 1 sweep, TLARS, TAFS-0.3, and TOMP track each other closely: FDR sits around 0.45–0.55 at SNR 0.01 with near-zero TPR, and by SNR >= 1 all three fall to FDR ~0.03–0.09 with TPR climbing to 1.00 at SNR 2–5. This matches the expectation that the screening decision depends mainly on `Phi_j` rather than on solver-specific path details — differences are within Monte Carlo noise at 20 runs. The Estimated FDR tracks the realized FDR and sits somewhat above it at high SNR (~0.09–0.10).
- **The 6-series part is where a solver interaction appears — via the bootstrap step, not the path algorithm.** The Ordinary rows in Part 2 reproduce Part 1. Under Bootstrap-CI, TLARS and TAFS-0.3 become markedly more conservative: lower FDR but a clear TPR cost (TLARS TPR 0.32 at SNR 1 and 0.90 at SNR 5; TAFS-0.3 0.72 and 0.75), whereas `TOMP (Bootstrap)` behaves almost like its Ordinary self (TPR 0.89 at SNR 1, 1.00 at SNR 2–5, FDR ~0.06). So the Ordinary-vs-Bootstrap gap is *not* uniform across solvers — TOMP interacts differently with the bootstrap thresholding step than the LARS-style backends do.
- `TAFS-0.3`'s extra regularization (`rho_afs = 0.3`) is the most likely place a solver-specific difference could surface, since it alters the T-loop's adaptive forward-stagewise behavior rather than just the path algorithm; at 20 runs its Ordinary curve stays within noise of TLARS, so any effect here is small.

**Last updated**: 2026-07-08
