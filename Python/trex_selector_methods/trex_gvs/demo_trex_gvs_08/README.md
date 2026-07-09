# Demo 08: T-Rex+GVS Block Benchmark — HAC-Discovered vs Oracle Groups (Python)

## Purpose

Benchmark how much accurate **group discovery** (via HAC clustering) costs
compared to knowing the true groups in advance ("oracle" groups), across three
within-block truth patterns. Structurally the most different demo in the suite:
instead of an SNR/rho sweep for a single method set, it runs a **4-method x
3-scenario** comparison over a `(rho x snr)` grid.

Python port of `cpp/.../demo_trex_gvs_08_mc_sim_block_bench`. The R counterpart is
`R/trex_selector_methods/trex_gvs/demo_trex_gvs_08/`.

Indices are **0-based** throughout (Python convention); the R counterpart is
1-based. `prior_groups` are contiguous 0-based block labels `0..G-1`, and
`active_blocks` are 0-based block IDs.

## Data-generating process

`dgp_block_equicorr` (from the shared `trex_gvs_dgps` module) — `G = 100`
equal-size blocks of `block_size = 5` (`p = 500`), with `n_active_blocks = 10`
chosen at random; each active block gets a random-sign coefficient of magnitude
`b = 3.0`. Columns follow `X_j = sqrt(rho) * Z_block(j) + sqrt(1 - rho) * eps`,
Gaussian, unit column variance. `n = 200`. The response uses the biased
(n-denominator) signal variance, matching the cpp benchmark DGP. Three truth
scenarios control which variables inside an active block are active:

- **INDIVIDUAL** — 1 active var/block (`s = 10`).
- **REPRESENTATIVE** — 2-3 active vars/block (`s ~ 20`-30).
- **WHOLE** — all `block_size` (5) vars/block (`s = 50`).

Grid: `rho_grid = {0.5, 0.9}`, `snr_grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`,
`num_MC = 500` per cell. Note `corr_max = 0.5` here — much lower than the `0.98`
used in demos 01-07, reflecting the deliberately weaker within-block correlation.
The per-trial data seed is staggered per grid cell
(`cell_base_seed = base_seed + (ir * n_snr + isnr) * num_MC`, `base_seed = 2026`),
so each `(rho, snr)` cell draws from a distinct seed band; the same `trial_seed` is
reused as the selector `cv_seed`.

## Methods compared

Four variants = EN/IEN x HAC-discovered/oracle groups. EN variants use
`en_solver = "TENET_AUG"`; IEN variants use `TENET`:

- **M1 (EN-C)** — `gvs_type = "EN"`, HAC-clustered groups (discovered from the
  data, `corr_max = 0.5`, single linkage).
- **M2 (EN-O)** — `gvs_type = "EN"`, oracle groups (the DGP's true block labels,
  supplied via the GVS control's `prior_groups` field as contiguous 0-based
  labels).
- **M3 (IE-C)** — `gvs_type = "IEN"`, HAC-clustered groups.
- **M4 (IE-O)** — `gvs_type = "IEN"`, oracle groups.

Shared controls (from `CFG`): `tFDR = 0.1`, `K = 20`, single linkage,
`lambda2_method = "CV_1SE_CCD"`. The central comparisons are **M1-vs-M2** and
**M3-vs-M4**: how much accuracy automatic HAC group discovery costs relative to
knowing the true groups.

## Parts

- **Part 1** — Single-run demo: all four variants on one INDIVIDUAL-scenario draw
  (`rho = 0.9`, `snr = 2.0`), printing a small metrics table with TPR, FDR,
  block_FDR, blk_hit, and T_stop.
- **Part 2** — Full benchmark over the three truth scenarios INDIVIDUAL /
  REPRESENTATIVE / WHOLE, each `(rho, snr)` cell averaged over `num_MC = 500`
  trials and printed as one grid table per `rho` value.

## Metrics

Averaged across the `num_MC` trials of each cell:

- **coord TPR / FDR** — coordinate-level (per-variable) true-positive and
  false-discovery proportions.
- **block_FDR** — block-level FDP = (null blocks hit) / max(total blocks hit, 1).
- **scenario-specific block metric** — `blk_hit` for INDIVIDUAL / REPRESENTATIVE
  (fraction of active blocks with `>= 1` selected variable); `full_blk` for WHOLE
  (fraction of active blocks fully selected).
- **T_stop** — the selector's stopping dummy count.
- **purity** — block-purity rate: the fraction of the `G` true blocks whose members
  all share one discovered-cluster label. `1.0` for oracle runs by construction;
  `< 1.0` for HAC runs indicates the clustering split or merged true blocks.
- **M_found** — number of clusters the selector formed.

### Grouping diagnostics — `purity` and `M_found`

Both grouping diagnostics are read from the `GVSSelectionResult` returned by
`select()`, which exposes `max_clusters` (→ `M_found`) and `groups_vec` (the
per-variable cluster labels). `purity` is computed from `groups_vec` via the same
block-purity rule as the cpp benchmark (`_block_purity_rate` in the demo), for HAC
**and** oracle methods alike; `M_found` is `max_clusters` for every method.

This actually surfaces *more* than the R demo: the R6 API exposes only
`$max_clusters`, so the R demo reports `purity` as `1.0` (oracle) / `-` (HAC)
rather than the real per-block rate. The Python `GVSSelectionResult` carries the
per-variable labels, so the genuine HAC purity is available here — e.g. at
`corr_max = 0.5` it falls toward `0` as `rho` drops (blocks shatter) and rises to
`1.0` at high `rho` (blocks recovered cleanly).

## Running

Only **Part 1 runs by default**; the full benchmark (Part 2) is opt-in via the
`RUN_PART_2` flag at the top of the file. (This differs from the R suite, where all
parts default on.)

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.py      # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.py 8    # 8 worker cores
```

The optional positional argument sets the number of Monte Carlo worker cores
(default `min(6, os.cpu_count())`). The file inserts both its own directory and the
parent suite directory onto `sys.path`, so the shared `trex_gvs_dgps` and
`gvs_sim_common` modules resolve and remain re-importable by the `spawn`-launched
workers (the top-level `_block_cell_worker` is imported by each worker) — run the
demo as a file, not via `python -`. Output is console-only (no result files, no
Plotly heatmap is rendered).

## Interpretation

- Central comparison is M1-vs-M2 (and M3-vs-M4): if `corr_max` is well matched to
  the true within-block `rho`, HAC-clustered groups should track oracle-group
  performance closely; a persistent gap indicates the clustering threshold is
  splitting or merging true blocks.
- INDIVIDUAL is expected to be hardest (only 1 of 5 correlated variables per block
  is truly active, so a discovered group can "win" on a correlated null neighbor);
  WHOLE should be easiest.
- Comparing `rho = 0.5` vs `rho = 0.9` shows how within-block correlation strength
  helps group discovery but can hurt localizing signal within a block in
  INDIVIDUAL / REPRESENTATIVE.
- EN vs IEN at matched group-source shows whether the iterative elastic net changes
  the HAC-vs-oracle gap itself or only the overall detection level.
- `purity` and `M_found` directly diagnose HAC grouping quality: for the HAC methods
  (M1, M3), purity near `1.0` with `M_found ≈ G` means the clustering recovered the
  true blocks; a low purity with a large `M_found` means it shattered them. Read
  these alongside `block_FDR` and the block-hit metrics.

**Last updated**: 2026-07-08
