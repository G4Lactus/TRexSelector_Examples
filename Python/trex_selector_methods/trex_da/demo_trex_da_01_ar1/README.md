# Demo 01: T-Rex+DA+AR1 on AR(1) Toeplitz Data (Python)

## Purpose

The foundational dependency-aware scenario: recover a sparse support under AR(1)
Toeplitz column correlation using the `AR1` DA method. A single-run illustration is
paired with Monte Carlo sweeps over SNR and rho, plus a 2D gap-vs-rho study that
probes exactly when the AR1 correction window starts suppressing true signals.

Python port of `R/trex_selector_methods/trex_da/demo_trex_da_01_ar1` and
`cpp/.../demo_trex_da_01_mc_sim_ar1`.

## Data-generating process

`dgp_ar1_snr` builds a Toeplitz design where `Sigma[j,k] = rho^|j-k|` via the
recursion `X[i,j] = rho * X[i,j-1] + sqrt(1 - rho^2) * eta[i,j]`, then calibrates
the noise to a target SNR. Single-run config: `n=150, p=500, s=5, amplitude=3.0,
rho=0.7, snr=2.0, tFDR=0.2, K=20, seed=2026`. Monte Carlo config: `n=300, p=1000,
s=10, amplitude=3.0, rho=0.7, snr=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026,
max_gap=20`.

## DA method

`build_da_selector("AR1", ...)` constructs a `trex_selector_neo.TRexDASelector`
with a `TRexDAControlParameter` (`method = DAMethod.AR1`, `rho_thr_DA = 0.02`,
`cor_coef` left auto) and a `TRexControlParameter` (`K = 20`,
`solver_type = SolverTypeForTRex.TLARS`).

## What it runs

Parts are toggled by the `RUN_PART_*` flags at the top of the file. Default: Part 1
only (`RUN_PART_1 = True`; Parts 2-4 `False`).

- Part 1: single-run demo on AR(1) data; prints true support, T_stop/L, selection, TPP/FDP.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` at rho=0.7, CappedSpread(max_gap=20) vs Random (per-trial) support.
- Part 3: rho sweep `{0.0, 0.1, ..., 0.9}` at SNR=2.0, CappedSpread vs Random support.
- Part 4: 2D gap x rho sweep; `gap_grid = {100, 50, 20, 15, 10, 5, 1}`, `rho_grid = {0.0..0.9}`, with the kappa boundary `kappa = ceil(log(0.02)/log(rho))`.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Enable the
heavier Monte Carlo parts by flipping the `RUN_PART_*` flags. Output is console-only;
the folder carries its own `simulation_results/` directory for optionally recorded
summaries. Monte Carlo parts run in parallel via Python `multiprocessing` using the
`spawn` start method (up to `min(6, os.cpu_count())` workers).

## Interpretation

- The AR1 correction substantially reduces FDR relative to the classical (no-DA) T-Rex on the same correlated data, at a real TPR cost as SNR grows since some correlated neighbors of true positives are actively excluded.
- The 2D gap x rho study is the key diagnostic: once active variables sit closer together than the kappa boundary implied by rho, they fall inside each other's correction windows and TPP collapses. This is a structural property, not a bug.

**Last updated**: 2026-07-08
