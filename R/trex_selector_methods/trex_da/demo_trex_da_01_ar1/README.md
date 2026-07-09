# Demo 01: T-Rex+DA+AR1 on AR(1) Toeplitz Data

## Purpose

The foundational dependency-aware scenario: recover a sparse support under AR(1)
Toeplitz column correlation using the `AR1` DA method. The demo pairs a single-run
illustration with Monte Carlo sweeps over SNR and rho, plus a 2D gap-vs-rho study
that probes exactly when the correction window starts suppressing true signals.

## Data-generating process

`dgp_ar1_snr` builds a Toeplitz design where `Sigma[j,k] = rho^|j-k|` via the
recursion `X[i,j] = rho * X[i,j-1] + sqrt(1 - rho^2) * eta[i,j]`, then calibrates
noise to a target SNR. Single-run config: `n=150, p=500, s=5, amplitude=3.0,
rho=0.7, SNR=2.0`. Monte Carlo config: `n=300, p=1000, s=10, amplitude=3.0,
rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "AR1",
rho_thr_DA = 0.02), control = trex_control(solver = "TLARS", K = 20))`. The
AR1 method estimates the correlation coefficient automatically (`cor_coef = NA`).

## What it runs

- Part 1: single-run demo on AR(1) data, prints selection/TPP/FDP and a correlation heatmap.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`, CappedSpread(max_gap=20) vs Random support.
- Part 3: rho sweep `seq(0.0, 0.9, by = 0.1)` at SNR=2.0, CappedSpread vs Random support.
- Part 4: 2D gap x rho sweep; `gap_grid = {100, 50, 20, 15, 10, 5, 1}`, `rho_grid = {0.0..0.9}`, with the kappa-boundary `kappa = ceiling(log(0.02)/log(rho))`.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- The AR1 correction substantially reduces FDR relative to the classical (no-DA) T-Rex on the same correlated data, at a real TPR cost as SNR grows since some correlated neighbors of true positives are actively excluded.
- The 2D gap x rho study is the key diagnostic: once active variables sit closer together than the kappa-boundary implied by rho, they fall inside each other's correction windows and TPP collapses. This is a structural property, not a bug.
- Comparing CappedSpread vs Random support shows whether deterministic even spacing behaves differently from random placement at the same nominal spacing.

**Last updated**: 2026-07-08
