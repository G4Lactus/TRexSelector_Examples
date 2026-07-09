# Demo 03b: T-Rex+DA+NN on AR(1) Data — Method-Mismatch Study (Python)

## Purpose

A deliberate misspecification stress test: apply the `NN` (banded / MA-band)
correction to AR(1) Toeplitz data — a geometrically-decaying, not banded,
correlation structure — to probe how the finite-range NN correction behaves when the
true dependence never truly vanishes at any distance. Companion to Demo 03 (`NN`
correction on correctly-specified banded data).

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_04_nn_02_ar` and
`cpp/.../demo_trex_da_03b_mc_sim_nn_ar`.

## Data-generating process

`dgp_ar1_snr` (the same AR(1) Toeplitz DGP as Demo 01): `X[i,j] = rho * X[i,j-1] +
sqrt(1 - rho^2) * eta[i,j]`. Single-run config: `n=150, p=500, s=5, rho=0.7,
snr=2.0, amplitude=3.0, tFDR=0.2, K=20, seed=2026`. Monte Carlo config: `n=300,
p=1000, s=10, amplitude=3.0, rho=0.7, snr=2.0, tFDR=0.2, K=20, num_MC=200,
seed=2026, max_gap=20`. Sweeps use CappedSpread(max_gap=20) support.

## DA method

`build_da_selector("NN", ...)` constructs a `trex_selector_neo.TRexDASelector` with
a `TRexDAControlParameter` (`method = DAMethod.NN`, no extra DA arguments) and a
`TRexControlParameter` (`K = 20`, `solver_type = SolverTypeForTRex.TLARS`). The NN
method is intentionally applied to AR(1) data here rather than the matched `AR1`
method of Demo 01.

## What it runs

Parts are toggled by the `RUN_PART_*` flags. Default: Part 4 only
(`RUN_PART_4 = True`; Parts 1-3 `False`), matching the R demo.

- Part 1: single-run demo (NN correction on AR(1) data, rho=0.7, SNR=2.0).
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` at rho=0.7, CappedSpread support.
- Part 3: rho sweep `{0.1, ..., 0.9}` at SNR=2.0, CappedSpread support.
- Part 4: 2D SNR x rho sweep; SNR `{2.0, 5.0}` x rho `{0.1, ..., 0.9}`, CappedSpread support.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_03b_nn_ar/demo_trex_da_03b_nn_ar.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Toggle
parts via the `RUN_PART_*` flags. Output is console-only; the folder carries its own
`simulation_results/` directory for optionally recorded summaries. Monte Carlo parts
run in parallel via Python `multiprocessing` using the `spawn` start method (up to
`min(6, os.cpu_count())` workers).

## Interpretation

- This tests robustness to model misspecification: NN assumes correlation vanishes beyond kappa columns, while AR(1) correlation decays but never reaches zero. Expect NN to recover some of the FDR-control benefit of the matched `AR1` correction (Demo 01) but likely less completely, as residual correlation beyond the NN window leaks through as extra false discoveries.
- Read side-by-side with Demo 01 (matched AR1) and Demo 03 (matched NN) to triangulate how much of DA-TRex's benefit comes from the exactly-right dependency model versus an approximately reasonable one.

**Last updated**: 2026-07-08
