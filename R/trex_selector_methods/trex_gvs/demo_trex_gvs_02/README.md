# Demo 02: T-Rex+GVS on Scattered Blocks

## Purpose

Test whether spatially **scattering** the members of a correlated group across the design matrix — rather than keeping them contiguous — affects T-Rex+GVS group recovery. This isolates the effect of column layout from correlation strength: the group structure is identical to Demo 01, only the column positions are permuted.

## Data-generating process

`dgp_scattered_grouped_snr` — `n_groups = 3` equicorrelated latent-factor groups of `group_size = 50` variables each (all active, `beta = 3`, `s = 150`), whose members are assigned by a random column permutation instead of contiguous blocks. Background columns are i.i.d. `N(0,1)`. Each grouped column is `X_j = Z_g + sd_x * xi`, with within-group correlation `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`.

## Methods compared

- **EN** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `GVS_type = "IEN"`, `en_solver = "TENET"`

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: correlation heatmap (scattered off-diagonal structure) + one EN / EN+AUG / IEN run.
- **Part 2** — SNR sweep at fixed `sd_x = sqrt(0.01)` (`rho ~ 0.99`); `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`; `rho in {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho).

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_02/demo_trex_gvs_02.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_02/demo_trex_gvs_02.R 8    # 8 worker cores
```

Output is console-only (no result files are written); the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_gvs_dgps.R` and `../../simulation_utils.R`.

## Interpretation

- Because HAC-based grouping clusters columns by **correlation**, not by column adjacency, scattering group members is expected to leave FDP/TPP patterns very close to Demo 01's contiguous-block results.
- Best read side-by-side with Demo 01: similar results confirm that group discovery depends on correlation structure rather than spatial layout; a noticeable divergence would point to a layout-sensitive bug or clustering limitation.
- As in Demo 01, expect IEN to remain conservative (very low FDP and TPP) on this all-active-groups design.

**Last updated**: 2026-07-08
