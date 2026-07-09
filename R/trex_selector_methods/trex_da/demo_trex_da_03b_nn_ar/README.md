# Demo 03b: T-Rex+DA+NN on AR(1) Data (Method-Mismatch Study)

## Purpose

A deliberate misspecification stress test: apply the `NN` (banded/nearest-neighbor)
DA correction to data generated from an AR(1) process, a genuinely different
(geometrically-decaying, not banded) correlation structure. It asks whether the NN
correction, which assumes a finite-range dependency, can still usefully approximate
a structure whose correlation never truly vanishes. Companion to Demo 03 (NN on
correctly-specified NN data).

## Data-generating process

`dgp_ar1_snr` builds the same AR(1) Toeplitz design as Demo 01
(`X[i,j] = rho * X[i,j-1] + sqrt(1 - rho^2) * eta[i,j]`), but here the selector is
configured with the NN method instead of AR1. Single-run config: `n=150, p=500,
s=5, amplitude=3.0, rho=0.7, SNR=2.0`. Monte Carlo config: `n=300, p=1000, s=10,
amplitude=3.0, rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "NN"),
control = trex_control(solver = "TLARS", K = 20))` — the NN correction applied to
AR(1) data (`trex_method = "trex+DA+NN"` in the MC runner).

## What it runs

- Part 1: single-run demo, NN correction on AR(1) data, prints selection/TPP/FDP and a correlation heatmap.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`, CappedSpread(max_gap=20) vs Random support.
- Part 3: rho sweep `seq(0.1, 0.9, by = 0.1)` at SNR=2.0, CappedSpread vs Random support.
- Part 4: 2D SNR x rho sweep; `SNR in {2.0, 5.0}`, `rho in seq(0.1, 0.9, by = 0.1)`, Random support (the narrow operating region of DA+NN under AR(1)).

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03b_nn_ar/demo_trex_da_03b_nn_ar.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03b_nn_ar/demo_trex_da_03b_nn_ar.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- This directly tests robustness to model misspecification: NN assumes correlation vanishes beyond kappa columns, while AR(1) correlation decays but never truly reaches zero, so some residual correlation should leak through as extra false discoveries relative to Demo 01's matched AR1-on-AR1 results.
- Read side-by-side with Demo 01 (same DGP, correctly-specified AR1) and Demo 03 (correctly-specified NN, different DGP) to gauge how much of the benefit comes from having exactly the right dependency model vs. an approximately reasonable one.
- If DA-NN performs nearly as well as DA-AR1 here, the method is reasonably robust to this mismatch; a large gap would indicate sensitivity to matching the true dependency structure.

**Last updated**: 2026-07-08
