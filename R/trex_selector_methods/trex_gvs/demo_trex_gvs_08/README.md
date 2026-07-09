# Demo 08: T-Rex+GVS Block Benchmark — HAC-Discovered vs Oracle Groups

## Purpose

Benchmark how much accurate **group discovery** (via HAC clustering) costs compared to knowing the true groups in advance ("oracle" groups), across three within-block truth patterns. Structurally the most different demo in the suite: instead of an SNR/rho sweep for a single method set, it runs a **4-method x 3-scenario** comparison over a `(rho x snr)` grid.

## Data-generating process

`dgp_block_equicorr` — `G = 100` equal-size blocks of `block_size = 5` (`p = 500`), with `n_active_blocks = 10` chosen at random; each active block gets a random sign coefficient of magnitude `b = 3.0`. Columns follow `X_j = sqrt(rho) * Z_g + sqrt(1 - rho) * eps`, Gaussian. `n = 200`. Three truth scenarios control which variables inside an active block are active: **INDIVIDUAL** (1/block, `s = 10`), **REPRESENTATIVE** (2–3/block, `s ~ 20`–30), **WHOLE** (all 5/block, `s = 50`).

Grid: `rho_grid = {0.5, 0.9}`, `snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`, `num_MC = 500` per cell. Note `corr_max = 0.5` here (vs 0.98 in demos 01–07). `cv_seed` is staggered per grid cell (`cell_base_seed = base_seed + cell_index * num_MC`), so each `(rho, snr)` cell draws from a distinct seed band (reused across the three scenarios).

## Methods compared

Four variants = EN/IEN x HAC-discovered/oracle groups:

- **M1 (EN-C)** — `GVS_type = "EN"`, HAC-clustered groups. EN variants use `en_solver = "TENET_AUG"`.
- **M2 (EN-O)** — `GVS_type = "EN"`, oracle groups (the DGP's true block labels, passed via `trex_gvs_control(groups = ...)` as contiguous 1..G IDs — GVS uses `groups`, not the DA-only `prior_groups`).
- **M3 (IE-C)** — `GVS_type = "IEN"`, HAC-clustered groups.
- **M4 (IE-O)** — `GVS_type = "IEN"`, oracle groups.

Shared controls: `tFDR = 0.1`, `K = 20`, `hc_linkage = "Single"`, `lambda2_method = "CV_1SE_CCD"`. The central comparisons are M1-vs-M2 and M3-vs-M4.

## Parts

- **Part 1** — Single-run demo: all four variants on one INDIVIDUAL-scenario draw (`rho = 0.9`, `snr = 2.0`), reporting TPR/FDR, block_FDR, blk_hit, T_stop, M_found.
- **Part 2** — Full benchmark over the three truth scenarios INDIVIDUAL / REPRESENTATIVE / WHOLE, each printed as one grid table per rho value.

Metrics: coordinate TPR/FDR, `block_FDR`, a scenario-specific block metric (`blk_hit` for INDIVIDUAL/REPRESENTATIVE, `full_blk` for WHOLE), `T_stop`, and `M_found`.

**Note on `purity` / `M_found`**: the R6 `TRexSelectorNeo` API does not expose the selector's internal per-variable cluster labels, so the C++ "purity" diagnostic is reported as `1.0` for the oracle methods (M2/M4, true by construction) and `-` (not available) for the HAC methods (M1/M3). Use `M_found` (number of clusters formed, read from `$max_clusters`) as the R-side grouping diagnostic instead.

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.R 8    # 8 worker cores
```

Output is console-only (no result files are written); the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_gvs_dgps.R`.

## Interpretation

- Central comparison is M1-vs-M2 (and M3-vs-M4): if `corr_max` is well matched to the true within-block `rho`, HAC-clustered groups should track oracle-group performance closely; a persistent gap indicates the clustering threshold is splitting or merging true blocks.
- INDIVIDUAL is expected to be hardest (only 1 of 5 correlated variables per block is truly active, so a discovered group can "win" on a correlated null neighbor); WHOLE should be easiest.
- Comparing `rho = 0.5` vs `rho = 0.9` shows how within-block correlation strength helps group discovery but can hurt localizing signal within a block in INDIVIDUAL/REPRESENTATIVE.
- EN vs IEN at matched group-source shows whether the iterative elastic net changes the HAC-vs-oracle gap itself or only the overall detection level.

**Last updated**: 2026-07-08
