# Demo 02: T-Rex+GVS on Scattered Blocks (Python)

## Purpose

Test whether spatially **scattering** the members of a correlated group across the design matrix — rather than keeping them contiguous — affects T-Rex+GVS group recovery. This isolates the effect of column layout from correlation strength: the group structure is identical to Demo 01, only the column positions are permuted.

Python port of `cpp/.../demo_trex_gvs_02_mc_sim_scattered_blocks`. The R counterpart is `R/trex_selector_methods/trex_gvs/demo_trex_gvs_02/`. Column indices here are **0-based** (Python convention); the R/cpp counterparts are 1-based.

## Data-generating process

`dgp_scattered_grouped_snr` (from `trex_gvs_dgps.py`) — `n_groups = 3` equicorrelated latent-factor groups of `group_size = 50` variables each (all active, `beta = 3`, `s = 150`), whose members are assigned by a random column permutation instead of contiguous blocks. Background columns are i.i.d. `N(0,1)`. Each grouped column is `X_j = Z_g + sd_x * xi`, with within-group correlation `rho = 1 / (1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls (`CFG`): `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: one EN / EN+AUG / IEN run, printing each result block. No correlation heatmap is rendered (unlike the R Part 1).
- **Part 2** — SNR sweep at fixed `sd_x = sqrt(0.01)` (`rho ~ 0.99`); `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`; `rho in {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho).

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_02/demo_trex_gvs_02.py      # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_02/demo_trex_gvs_02.py 8    # 8 worker cores
```

Only **Part 1** (the single run) runs by default. The heavy Monte Carlo parts (200 MC per grid point) are **opt-in**: flip the corresponding `RUN_PART_2 … RUN_PART_4` flags at the top of the file to `True`. This differs from the R suite, where all parts default on.

The demo inserts both its own directory and the parent suite directory onto `sys.path`, so the shared `trex_gvs_dgps` / `gvs_sim_common` modules resolve from the nested location and are re-importable by the `spawn`-launched Monte Carlo workers. Run it **as a file** (not via `python -` / piped stdin), or the workers cannot re-import the modules. The optional first CLI argument overrides the worker-core count (default `min(6, os.cpu_count())`). Output is console-only — no result files are written.

## Interpretation

- Because HAC-based grouping clusters columns by **correlation**, not by column adjacency, scattering group members is expected to leave FDP/TPP patterns very close to Demo 01's contiguous-block results.
- Best read side-by-side with Demo 01: similar results confirm that group discovery depends on correlation structure rather than spatial layout; a noticeable divergence would point to a layout-sensitive bug or clustering limitation.
- As in Demo 01, expect IEN to remain conservative (very low FDP and TPP) on this all-active-groups design.

**Last updated**: 2026-07-08
