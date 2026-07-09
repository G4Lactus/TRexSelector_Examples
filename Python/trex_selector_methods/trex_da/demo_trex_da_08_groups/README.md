# Demo 08: T-Rex+DA on Multi-Level Prior-Groups Data (Python)

## Purpose

Test DA-TRex on the structurally richest dependency model in this suite: a
three-level hierarchical latent-factor structure with a non-exchangeable Toeplitz
leaf layer, using explicit prior group labels rather than an automatically-discovered
dendrogram cut. This isolates the difficulty of the multi-level dependency-aware
correction itself from any clustering-quality artifact.

Python port of `cpp/.../demo_trex_da_08_mc_sim_groups`. The R counterpart is
`R/trex_selector_methods/trex_da/demo_trex_da_08_groups/`: the TRexSelectorNeo R
binding's `trex_da_control()` gained `prior_groups` / `rho_grid_labels` in the
2026-07 upstream fix, so the prior-groups method is now drivable from R as well.

## Data-generating process

`dgp_groups_toeplitz_leaf` draws from a multi-level latent-factor model where coarser
levels use shared factors and the finest level uses a Toeplitz-correlated leaf block
(`Cov(u_r, u_s) = phi^|r-s|`). Three nested grouping levels: groups of 10, 50, and
250 (fine to coarse), built by `_group_structure(p)` as `[j//10, j//50, j//250]`.
Config: `n=300, p=1000, s=10, amplitude=1.0, tFDR=0.2, K=20, num_MC=200, seed=2026`;
cumulative correlation levels `rho_levels = {0.55, 0.25, 0.10}`, Toeplitz leaf decay
`phi_leaf = 0.50`, Random support redrawn per trial.

## DA method

The prior-groups path: `build_da_selector("PRIOR_GROUPS", ..., prior_groups=groups,
rho_grid_labels=[0.55, 0.25, 0.10])` sets `TRexDAControlParameter.method =
DAMethod.PRIOR_GROUPS`, `prior_groups` (the 3-level hierarchy), and `rho_grid_labels`.
Because the hierarchy is supplied directly (non-empty `prior_groups`), it is used as-is
rather than being discovered by HAC/BT clustering. Solver is
`SolverTypeForTRex.TLARS`, `K = 20`.

## What it runs

Part 1 active by default (`RUN_PART_1 = True`).

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` (n=300, p=1000, s=10, Random support redrawn per trial), 200 MC per point.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_da/demo_trex_da_08_groups/demo_trex_da_08_groups.py
```

The file inserts both its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` and `da_sim_common` modules import directly. Output is
console-only; the folder carries its own `simulation_results/` directory for
optionally recorded summaries. Monte Carlo runs in parallel via Python
`multiprocessing` using the `spawn` start method (up to `min(6, os.cpu_count())`
workers).

## Interpretation

- The three-level nested dependency structure is genuinely harder for DA-TRex to control than the single-level structures in Demos 01-07, so realized FDR can overshoot the tFDR=0.2 target more than elsewhere in this suite.
- Because prior groups are supplied directly rather than HAC-discovered, any residual FDR inflation reflects genuine difficulty in the underlying multi-level correction, not a clustering-quality confound — which differentiates Demo 08 from the BT-based Demos 04-07. The demo currently exercises only the SNR sweep.

**Last updated**: 2026-07-08
