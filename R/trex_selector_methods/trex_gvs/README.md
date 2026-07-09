# T-Rex+GVS: Grouped Variable Selection — R Demonstration Suite

> **Status**: migrated to the **TRexSelectorNeo** R6 API
> (`TRexGVSSelector$new(..., gvs_control = trex_gvs_control(...))$select()`,
> with `compute_fdp()` / `compute_tpp()`). No longer depends on the CRAN
> `TRexSelector` package.

## Purpose

R demos and Monte Carlo simulations for the **Grouped Variable Selection
T-Rex (T-Rex+GVS)** selector on designs where the active variables cluster
into groups. The folder layout mirrors the C++ suite one-to-one: one
subfolder per demo, plus a single consolidated DGP file.

```txt
trex_gvs/
  ├── README.md
  ├── trex_gvs_dgps.R                 <- all scenario DGPs (one file)
  ├── demo_trex_gvs_01/
  │   ├── demo_trex_gvs_01.R
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_gvs_02/
  │   └── ...
  └── demo_trex_gvs_08/
```

Each demo begins with a **Part 1 single-run** (data/structure diagnostics plus
a Plotly correlation heatmap in interactive sessions) and then runs Monte Carlo
sweeps that compare the three dummy-generation methods
**EN (TENET) / EN+AUG (TENET_AUG) / IEN**, sweeping SNR, a correlation knob,
and a 2-D SNR × correlation grid. Typical fixed settings: n = 200, p = 500,
tFDR = 0.1, K = 20, 200 MC runs (demo 08 uses 500).

Dependencies: `TRexSelectorNeo`, `plotly`, `parallel`. Each demo sources the
consolidated [trex_gvs_dgps.R](trex_gvs_dgps.R) (all scenario generators in one
file, mirroring the C++ `trex_gvs_dgps.hpp`) and the shared
[../support_generators.R](../support_generators.R) /
[../simulation_utils.R](../simulation_utils.R) (which provides the `.run_mc_*()`
EN/EN+AUG/IEN wrappers and the `.gvs_select()` selector adapter). Demo 08 is
self-contained (its block runner, metrics, and grid printer live in the demo
file). Output is console-only; each demo's `simulation_results/` folder holds
any manually-recorded result summaries (see demo 01).

---

## Demos

| # | Scenario | Sweeps | DGP (in `trex_gvs_dgps.R`) | cpp counterpart |
|---|---|---|---|---|
| [01](demo_trex_gvs_01/) | Hastie: 3 equicorrelated active groups of 50 (all active, s = 150) | SNR, rho, 2-D, **lambda2_method**, **hc_linkage** | `dgp_hastie_snr` | `demo_trex_gvs_01_mc_sim_hastie_en_blocks` |
| [02](demo_trex_gvs_02/) | Scattered-grouped: same 3 groups, members permuted across columns | SNR, rho, 2-D | `dgp_scattered_grouped_snr` | `demo_trex_gvs_02_mc_sim_scattered_blocks` |
| [03](demo_trex_gvs_03/) | Mixed blocks: 3 active blocks (20/50/80) + 1 trap (65), random order, noise gaps | SNR, rho, 2-D — **two presets** (default sd_x; fixed rho = 0.75) | `dgp_mixed_blocks_snr` | `demo_trex_gvs_03_mc_sim_mixed_blocks` |
| [04](demo_trex_gvs_04/) | Negative traps: active +Z/−Z group + two inactive +Z/−Z traps (s = 100) | SNR, rho, 2-D | `dgp_neg_traps_snr` | `demo_trex_gvs_04_mc_sim_neg_traps` |
| [05](demo_trex_gvs_05/) | Heavy-tailed blocks: demo-03 geometry with Student-t(3) factors/noise | SNR (fixed rho = 0.75), rho, 2-D | `dgp_t3_equi_blocks_snr` | `demo_trex_gvs_05_mc_sim_t3_blocks` |
| [06](demo_trex_gvs_06/) | AR(1) blocks: demo-03 geometry with within-block Cor = rho^\|j−k\| | SNR (fixed rho = 0.85), rho, 2-D | `dgp_ar1_blocks_snr` | `demo_trex_gvs_06_mc_sim_ar1_blocks` |
| [07](demo_trex_gvs_07/) | ARMA blocks: AR(1) / MA(3) / ARMA(2,1) active blocks + AR(1) trap | SNR (fixed ar_coef = 0.8), **ar_coef**, 2-D | `dgp_arma_blocks_snr` | `demo_trex_gvs_07_mc_sim_arma_blocks` |
| [08](demo_trex_gvs_08/) | Block benchmark: 100 blocks × 5, 10 active; HAC-discovered vs oracle groups | 4 methods (M1–M4) × 3 truth scenarios over a (rho × snr) grid | `dgp_block_equicorr` | `demo_trex_gvs_08_mc_sim_block_bench` |

Each demo folder has its own `README.md`. Every demo exposes `run_part_*` flags
at the top to enable/disable individual parts. Output is console-only, as in C++
for demo 08.

### Demo 08 notes (block benchmark)

Demo 08 compares four variants — **M1** = EN + HAC-discovered groups, **M2** =
EN + oracle groups, **M3** = IEN + HAC, **M4** = IEN + oracle — across three
within-block truth scenarios (**INDIVIDUAL**: 1 active var/block;
**REPRESENTATIVE**: 2–3/block; **WHOLE**: all 5/block). Oracle groups are
supplied via `trex_gvs_control(groups = ...)` with contiguous 1..G labels; HAC
variants discover groups from the data (`corr_max = 0.5`, single linkage). It
reports coordinate TPR/FDR, block-level FDR and a scenario-specific block metric
(`blk_hit` / `full_blk`), plus `T_stop` and `M_found`.

The C++ benchmark also reports a group-**purity** diagnostic computed from the
selector's internal per-variable cluster labels. The TRexSelectorNeo R6 API does
not expose those labels, so purity is reported as `1.0` for the oracle methods
(true by construction) and `-` for the HAC methods; **M_found** (the number of
clusters the selector formed, from `$max_clusters`) is the R-side grouping
diagnostic.

---

## DGPs

All scenario generators live in the single [trex_gvs_dgps.R](trex_gvs_dgps.R)
(mirroring the C++ `trex_gvs_dgps.hpp`): `dgp_hastie_snr`,
`dgp_scattered_grouped_snr`, `dgp_mixed_blocks_snr`, `dgp_neg_traps_snr`,
`dgp_t3_equi_blocks_snr`, `dgp_ar1_blocks_snr`, `dgp_arma_blocks_snr`
(+ the `.sim_arma()` helper), and `dgp_block_equicorr`. Each carries an
interactive `generate_*()` variant alongside the SNR-calibrated `*_snr()` one.
Earlier speculative, demo-less DGPs (factor-model, hapgen, neg-corr,
random/unequal blocks, the `dgp_gvs` collection) were removed in the
consolidation.

---

## Recorded results

Manually-recorded result summaries live in each demo's `simulation_results/`
folder, e.g.
[demo_trex_gvs_01/simulation_results/demo_trex_gvs_01_hastie_models.md](demo_trex_gvs_01/simulation_results/demo_trex_gvs_01_hastie_models.md).

---

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_01/demo_trex_gvs_01.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.R 8    # 8 worker cores
```

Scripts resolve their own directory, so they run from any working directory.
The optional first argument sets the number of parallel worker cores (default 6).
The `run_part_*` flags at the top of each file enable or disable individual
parts. These are full Monte Carlo simulations (200–500 trials on n = 200,
p = 500); reduce `num_MC` / grid sizes for a quick smoke run.

---

**Last updated**: 2026-07-08
