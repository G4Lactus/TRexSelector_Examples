# Demo 08: T-Rex+DA on Multi-Level Prior-Groups Data

## Purpose

Test DA-TRex on the structurally richest dependency model in this suite: a
three-level hierarchical latent-factor structure with a non-exchangeable Toeplitz
leaf layer, using explicit prior group labels rather than an automatically
discovered dendrogram cut. This isolates the difficulty of the multi-level
dependency-aware correction itself from any clustering-quality artifact.

## Data-generating process

`dgp_groups_toeplitz_leaf` draws from a multi-level latent-factor model where
coarser levels use shared factors and the finest level uses a Toeplitz-correlated
leaf block (`Cov(u_r, u_s) = phi^|r-s|`). Three nested grouping levels: groups of
10, 50, and 250 (fine to coarse). Config: `n=300, p=1000, s=10, amplitude=1.0,
tFDR=0.2, K=20, num_MC=200, seed=2026`; cumulative correlation levels
`rho_levels = {0.55, 0.25, 0.10}`, Toeplitz leaf decay `phi_leaf = 0.50`, Random
support.

## DA method

The prior-groups path: `trex_da_control(prior_groups = groups,
rho_grid_labels = c(0.55, 0.25, 0.10))` supplies the 3-level hierarchy explicitly.
The core routes on `prior_groups` being non-empty, so the `da_method` value is
irrelevant here (this is not the HAC/BT-discovery path). Solver is TLARS, `K = 20`.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` (n=300, p=1000, s=10, Random support redrawn per trial), 200 MC per point.

## Running

```bash
Rscript R/trex_selector_methods/trex_da/demo_trex_da_08_groups/demo_trex_da_08_groups.R      # full run
Rscript R/trex_selector_methods/trex_da/demo_trex_da_08_groups/demo_trex_da_08_groups.R 8    # 8 worker cores
```

Output is console-only; the script sources the consolidated `../trex_da_dgps.R` and
the shared `../../simulation_utils.R` + `../../support_generators.R`.

## Interpretation

- The three-level nested dependency structure is genuinely harder for DA-TRex to control than the single-level structures in Demos 01-07, so realized FDR can overshoot the tFDR=0.2 target more than elsewhere in this suite.
- Because prior groups are supplied directly rather than HAC-discovered, any residual FDR inflation reflects genuine difficulty in the underlying multi-level correction, not a clustering-quality confound — which differentiates Demo 08 from the BT-based Demos 04-07.
- The demo currently exercises only the SNR sweep; extending it with rho-level or group-size sweeps (analogous to Demos 04-07's Q/M sweeps) would help clarify which aspect of the nested structure drives the elevated FDR.

**Last updated**: 2026-07-08
