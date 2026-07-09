# Demo 01: T-Rex+GVS on Hastie Equicorrelated Blocks (Python)

## Purpose

Evaluate grouped-signal recovery under strong equicorrelation, using the classical three-group design of Zou & Hastie (2005). This is the canonical starting point for the T-Rex+GVS suite: every variable inside an active group is itself active, so the question is purely **group-level recovery**, not within-group sparse-support identification.

Python port of `cpp/.../demo_trex_gvs_01_mc_sim_hastie_en_blocks`. The R counterpart is `R/trex_selector_methods/trex_gvs/demo_trex_gvs_01/`. Column indices here are **0-based** (Python convention); the R/cpp counterparts are 1-based.

## Data-generating process

`dgp_hastie_snr` (from `trex_gvs_dgps.py`) — three equicorrelated latent-factor groups of 50 variables each (columns 0–49, 50–99, 100–149), all active with `beta = 3`, plus 350 i.i.d. background columns. Each column is `X_j = Z_g + sd_x * xi`, giving within-group correlation `rho = 1 / (1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`, `s = 150`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`; SNR is calibrated through `sigma_y`.

## Methods compared

- **EN** — `gvs_type = "EN"`, `en_solver = "TENET"`
- **EN+AUG** — `gvs_type = "EN"`, `en_solver = "TENET_AUG"`
- **IEN** — `gvs_type = "IEN"`, `en_solver = "TENET"`

Shared controls (`CFG`): `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — Single-run demo: one EN / EN+AUG / IEN run, printing each result block (selection counts, TPP/FDP). No correlation heatmap is rendered (unlike the R Part 1).
- **Part 2** — SNR sweep at fixed `sd_x = sqrt(0.01)` (`rho ~ 0.99`); `SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`.
- **Part 3** — rho sweep at fixed `SNR = 2.0`; `rho in {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).
- **Part 4** — 2-D SNR x rho grid (5 SNR x 6 rho).
- **Part 5** — `lambda2_method` comparison, `CV_1SE_SVD` vs `CV_1SE_CCD` (EN and IEN), fixed `SNR = 2.0`, `rho = 0.7`.
- **Part 6** — `hc_linkage` comparison, Single / Complete / Average / WPGMA (EN and IEN), same operating point, `lambda2 = CV_1SE_CCD`.

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_01/demo_trex_gvs_01.py      # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_01/demo_trex_gvs_01.py 8    # 8 worker cores
```

Only **Part 1** (the single run) runs by default. The heavy Monte Carlo parts (200 MC per grid point) are **opt-in**: flip the corresponding `RUN_PART_2 … RUN_PART_6` flags at the top of the file to `True`. This differs from the R suite, where all parts default on.

The demo inserts both its own directory and the parent suite directory onto `sys.path`, so the shared `trex_gvs_dgps` / `gvs_sim_common` modules resolve from the nested location and are re-importable by the `spawn`-launched Monte Carlo workers. Run it **as a file** (not via `python -` / piped stdin), or the workers cannot re-import the modules. The optional first CLI argument overrides the worker-core count (default `min(6, os.cpu_count())`). Output is console-only — no result files are written; the folder's `simulation_results/` directory holds a manually-recorded reference summary.

## Interpretation

- FDP should stay near the `tFDR = 0.1` target across the SNR grid for EN and EN+AUG, with TPP rising from weak recovery at low SNR toward near-complete recovery once `SNR >= 1`.
- IEN is expected to behave much more conservatively on this all-active design (very low FDP, slowly rising TPP); this pattern recurs across the other equicorrelated-block demos.
- The rho sweep should be fairly stable across correlation levels, since HAC-based grouping aggregates evidence across strongly correlated columns.
- Parts 5–6 check that the qualitative conclusions are not sensitive to the `lambda2` selection rule or the HAC linkage choice.

**Last updated**: 2026-07-08
