# Demo 03: T-Rex+GVS on Mixed Blocks (Python)

## Purpose

Evaluate T-Rex+GVS on a more realistic grouped design: four contiguous equicorrelated blocks of unequal size, placed in **random order** with **random white-noise gaps** between and around them, including one **inactive equicorrelated trap block**. A step up in realism from Demo 01's tidy contiguous Hastie blocks.

Python port of `cpp/.../demo_trex_gvs_03_mc_sim_mixed_blocks`. The R counterpart is `R/trex_selector_methods/trex_gvs/demo_trex_gvs_03/`. Column indices here are **0-based** (Python convention); the R/cpp counterparts are 1-based.

## Data-generating process

`dgp_mixed_blocks_snr` (from `trex_gvs_dgps.py`) — four equicorrelated blocks of sizes `{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one inactive trap block (size 65) that shares the active blocks' correlation structure. Blocks are shuffled into a random order and separated by random-width white-noise gaps within `p = 500`. Each block column is `X_j = Z_g + sd_x * xi`, `rho = 1 / (1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`. The realized block order and per-block 0-based column ranges are randomized per seed and printed by Part 1.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls (`CFG`): `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: block-layout diagnostics (realized block order, per-block 0-based column ranges, active/inactive flag) followed by one EN / EN+AUG / IEN run, printing each result block. No correlation heatmap is rendered (unlike the R Part 1).
- **Part 2** — **Mixed-Blocks** preset (default `sd_x = sqrt(0.01)`, `rho ~ 0.99`): SNR sweep + rho sweep (9 levels) + 2-D SNR x rho grid.
- **Part 3** — **Mixed-Blocks-Rho075** preset (fixed `rho = 0.75` for the SNR sweep): SNR sweep + a wider rho sweep (11 levels, adds 0.40 and 0.60) + 2-D grid whose rho axis uses 0.85 in place of 0.95.

Both presets share the SNR grid `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` and the 2-D SNR grid `{0.2, 0.5, 1.0, 2.0, 5.0}`. The two presets are defined in the `PRESETS` dict and run through the shared `_run_mixed_benchmark` helper.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_03/demo_trex_gvs_03.py      # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_03/demo_trex_gvs_03.py 8    # 8 worker cores
```

Only **Part 1** (the single run) runs by default. The heavy Monte Carlo presets (200 MC per grid point) are **opt-in**: flip `RUN_PART_2` (Mixed-Blocks) and/or `RUN_PART_3` (Mixed-Blocks-Rho075) at the top of the file to `True`. This differs from the R suite, where all parts default on.

The demo inserts both its own directory and the parent suite directory onto `sys.path`, so the shared `trex_gvs_dgps` / `gvs_sim_common` modules resolve from the nested location and are re-importable by the `spawn`-launched Monte Carlo workers. Run it **as a file** (not via `python -` / piped stdin), or the workers cannot re-import the modules. The optional first CLI argument overrides the worker-core count (default `min(6, os.cpu_count())`). Output is console-only — no result files are written.

## Interpretation

- FDP should stay controlled near `tFDR = 0.1` for EN / EN+AUG across both presets, since the inactive trap block should be excluded once its group is identified as null.
- TPP should climb with SNR; the random gaps and block order mainly test that the grouping logic is invariant to block placement rather than materially changing detection power.
- Comparing the default preset (`rho ~ 0.99`) against Rho075 (`rho = 0.75`) shows how performance degrades as within-block correlation decreases, all else equal.
- As in Demos 01–02, expect IEN to be substantially more conservative than the EN-based methods.

**Last updated**: 2026-07-08
