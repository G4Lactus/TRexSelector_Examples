# Demo 02: T-Rex+DA+EQUI on Equicorrelated Data

## Purpose

Cover the equicorrelation (compound-symmetry) dependency structure, where every
predictor is tied to the same shared latent factor, using the `EQUI` DA method.
Unlike AR(1), equicorrelation has no notion of column distance, so the correction
acts uniformly across all variables.

## Data-generating process

`dgp_bt_snr` is reused with `rho_within == rho_between`, which reduces the
binary-tree factor model to a pure equicorrelated design
`X[i,j] = sqrt(rho) * f[i] + sqrt(1 - rho) * eta[i,j]`. Single-run config:
`n=150, p=500, s=5, rho=0.25, SNR=2.0, n_blocks=5`. Monte Carlo config:
`n=300, p=1000, s=10, amplitude=3.0, rho_within=rho_between=0.25, tFDR=0.2,
K=20, num_MC=200, seed=2026, n_blocks=10`, Random support.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "EQUI"),
control = trex_control(solver = "TLARS", K = 20))`.

## What it runs

- Part 1: single-run EQUI demo (`rho=0.25`, SNR=2.0), prints the selected-variable count and a correlation heatmap.
- Part 2: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` at fixed `rho=0.25`, Random support, 200 MC per point.
- Part 3: 1D pure-equi rho sweep `seq(0.0, 0.9, by = 0.1)` at SNR=2.0 (sweeping rho with `rho_within == rho_between`), Random support.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_02_equi/demo_trex_da_02_equi.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_02_equi/demo_trex_da_02_equi.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- Because equicorrelation ties every column to one shared factor, expect DA-EQUI to behave qualitatively like the AR(1) case (Demo 01): reduced FDR vs. base T-Rex at some TPR cost, but the effect should be more uniform across variables.
- The rho sweep isolates how the strength of the shared factor drives selection difficulty: at rho=0 the design is effectively independent, and the correction has little to do.
- The SNR sweep shows the usual pattern of TPR rising toward a plateau while the DA correction keeps FDR near the tFDR=0.2 target.

**Last updated**: 2026-07-08
