# Demo 03: T-Rex+DA+NN on Banded (Nearest-Neighbor) Data

## Purpose

Study the `NN` DA method, which targets a banded MA(kappa) covariance where only
the kappa nearest-neighboring columns are correlated and correlation is exactly
zero beyond that band. This is the correctly-specified counterpart to Demo 03b
(NN correction applied to AR(1) data).

## Data-generating process

`dgp_nn_snr` builds `X[i,j] = sum_{l=0..kappa} theta_l * eta[i,j+l]` with
`theta_l ~ rho^l` normalized to unit variance, so `Sigma[j,k] = 0` once
`|j-k| > kappa`. Single-run config: `n=150, p=500, s=5, amplitude=3.0, kappa=5,
rho=0.7, SNR=2.0`. Monte Carlo config: `n=300, p=1000, s=10, amplitude=3.0,
kappa=3, rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "NN"),
control = trex_control(solver = "TLARS", K = 20))`.

## What it runs

- Part 1: single-run NN demo, prints selection/TPP/FDP and a correlation heatmap.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`, CappedSpread(max_gap=20) vs Random support.
- Part 3: rho sweep `seq(0.1, 0.9, by = 0.1)` at SNR=2.0, CappedSpread vs Random support.
- Part 4: kappa sweep `{1, 2, 3, 5, 7, 10, 15}` at rho=0.7, SNR=2.0, CappedSpread vs Random support.
- Part 5: 2D kappa x rho sweep; `kappa in {1,2,3,5,7,10,15}`, `rho in {0.1,0.3,0.5,0.7,0.8,0.9}` at SNR=2.0.
- Part 6: 2D SNR x rho sweep; `SNR in {0.1,0.5,1.0,2.0,5.0}`, `rho in {0.5,0.6,0.7,0.75,0.8,0.85,0.9}` at kappa=3.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03_nn/demo_trex_da_03_nn.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_03_nn/demo_trex_da_03_nn.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- Because NN correlation has a hard cutoff at distance kappa, DA-NN's correction is most effective when active variables are spaced further apart than kappa, analogous to Demo 01's gap-vs-kappa story but with a sharp rather than decaying profile.
- The kappa sweep (Part 4) and 2D kappa x rho sweep (Part 5) are the most informative views: as kappa grows the effective correlation neighborhood widens, and DA-NN needs progressively wider correction windows with a corresponding TPR cost.
- The SNR x rho sweep (Part 6) probes whether stronger signal can rescue the selector when the correlation geometry is unfavourable (high rho, densely packed actives).

**Last updated**: 2026-07-08
