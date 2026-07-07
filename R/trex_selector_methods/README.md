# T-Rex Selector Methods — R Demos

R demonstration suites for the T-Rex selector variants, mirroring
[cpp/trex_selector_methods/](../../cpp/trex_selector_methods/). Each variant
lives in its own sub-directory; Monte Carlo demos write results to a
source-side `simulation_results/` folder next to the demo file.

---

## Contents

| Folder | Status | What it covers |
|---|---|---|
| [trex/](trex/README.md) | NEW — TRexSelectorNeo R6 API | Classical T-Rex: toy-data run, single runs, MC studies (fixed/variable support), l-loop strategies, memory-mapped workflows, scalability benchmark (demos 00–07) |
| [trex_da/](trex_da/README.md) | NEW — TRexSelectorNeo R6 API | Dependency-Aware T-Rex (`TRexDASelector`) over AR(1), equicorrelated, nearest-neighbor, block-diagonal, and heavy-tailed designs (demos 01–08). The two NN demos skip cleanly until the R binding exposes `da_method="NN"` |
| [trex_gvs/](trex_gvs/README.md) | NEW — TRexSelectorNeo R6 API | Grouped Variable Selection T-Rex (`TRexGVSSelector`, EN vs IEN) over Hastie, scattered, mixed, negative-trap, equi-block, t-distributed, AR(1), and ARMA designs (demos 01–08 + simulation scripts) |
| [trex_spca/](trex_spca/README.md) | NEW — TRexSelectorNeo R6 API | T-Rex Sparse PCA MC evaluation (`TRexSPCASelector` / `TRexGVSSelector`) plus four C++/R cross-check probe scripts |
| [trex_screening/](trex_screening/README.md) | NEW — TRexSelectorNeo R6 API | Screen-TRex (`TRexScreeningSelector` / `TRexBiobankScreeningSelector`): in-memory and memory-mapped screening, correlated designs, biobank-scale, and solver comparison (6 demos) |

---

## Shared helpers

| File | Status | Description |
|---|---|---|
| [support_generators.R](support_generators.R) | Shared | Support-index construction helpers, sourced by demos across all variants: `make_support_capped_spread()` (deterministic evenly-spaced support with capped gap), `make_support_random()` (uniformly random support; internally uses `set.seed(seed + 500000L)`), `make_support_bt_one_per_block()` (one active predictor per BT block). All return 1-based R indices and mirror the C++ helpers. |
| [simulation_utils.R](simulation_utils.R) | NEW — migrated to the TRexSelectorNeo R6 API (`TRexDASelector` / `TRexGVSSelector` + `compute_fdp()` / `compute_tpp()`) | Shared MC engine for the trex_da and trex_gvs demos: `plot_cormat()` (Plotly heatmap), table printers, the generic parallel MC runner `.run_mc()`, and per-DGP wrappers (`.run_mc_ar1()`, `.run_mc_nn()`, `.run_mc_equi()`, `.run_mc_bt()`, the block-sweep wrappers, and the GVS wrappers `.run_mc_hastie()` … `.run_mc_hapgen_snr()` comparing EN vs IEN). |

The NEW-API trex demos additionally use their own utility file,
[trex/trex_sim_utils.R](trex/trex_sim_utils.R) (see the
[trex README](trex/README.md)).

---

## Running

All demos run from any working directory, e.g.

```bash
Rscript R/trex_selector_methods/trex/demo_trex_02_mc_sim_fixed_support.R
Rscript R/trex_selector_methods/trex_da/demo_trex_da_01_ar1.R
```

The da / gvs MC suites additionally require `plotly` (correlation heatmaps)
and use `parallel` for the MC loops; the spca demos additionally require
`elasticnet` and `glmnet`.

---

**Last updated**: 2026-07-08
