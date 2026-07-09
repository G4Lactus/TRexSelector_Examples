# Demo 03: T-Rex+GVS on Mixed Blocks

## Purpose

Evaluate T-Rex+GVS on a more realistic grouped design: four contiguous equicorrelated blocks of unequal size, placed in **random order** with **random white-noise gaps** between and around them, including one **inactive equicorrelated trap block**. A step up in realism from Demo 01's tidy contiguous Hastie blocks.

## Data-generating process

`dgp_mixed_blocks_snr` — four equicorrelated blocks of sizes `{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one inactive trap block that shares the active blocks' correlation structure. Blocks are shuffled into a random order and separated by random-width white-noise gaps within `p = 500`. Each block column is `X_j = Z_g + sd_x * xi`, `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`.

## Methods compared

- **EN** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `GVS_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: block-layout diagnostics (order, per-block column ranges) + one EN / EN+AUG / IEN run.
- **Part 2** — **Mixed-Blocks** preset (default `sd_x = sqrt(0.01)`, `rho ~ 0.99`): SNR sweep + rho sweep (9 levels) + 2-D SNR x rho grid.
- **Part 3** — **Mixed-Blocks-Rho075** preset (fixed `rho = 0.75` for the SNR sweep): SNR sweep + a wider rho sweep (11 levels, adds 0.40 and 0.60) + 2-D grid whose rho axis uses 0.85 in place of 0.95.

Both presets share the SNR grid `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` and 2-D SNR grid `{0.2, 0.5, 1.0, 2.0, 5.0}`.

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_03/demo_trex_gvs_03.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_03/demo_trex_gvs_03.R 8    # 8 worker cores
```

Output is console-only (no result files are written); the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_gvs_dgps.R` and `../../simulation_utils.R`.

## Interpretation

- FDP should stay controlled near `tFDR = 0.1` for EN / EN+AUG across both presets, since the inactive trap block should be excluded once its group is identified as null.
- TPP should climb with SNR; the random gaps and block order mainly test that the grouping logic is invariant to block placement rather than materially changing detection power.
- Comparing the default preset (`rho ~ 0.99`) against Rho075 (`rho = 0.75`) shows how performance degrades as within-block correlation decreases, all else equal.
- As in Demos 01–02, expect IEN to be substantially more conservative than the EN-based methods.

**Last updated**: 2026-07-08
