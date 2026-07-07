# TRexSelector_Examples — R Sources

## Overview

This directory contains the standalone R examples for the **TRexSelector** project.
It mirrors the sibling [cpp/](../cpp/README.md) examples tree and covers the same
library stack from the R side: memory-mapped matrices, data-preprocessing
building blocks (ml_methods), and T-Rex selector variants.

Two R packages appear in this tree:

- **TRexSelectorNeo** (v0.2.0) — the new Rcpp/R6 package from the main repo
  (`TRexSelector/R/TRexSelectorNeo`, requires R >= 4.3). Its R6 API is
  `TRexSelector$new(X, y, tFDR, control = trex_control(...))$select()`, plus
  helpers such as `compute_fdp()` / `compute_tpp()`, `mmap_matrix()` /
  `convert_to_memory_mapped()`, the scaler classes, and
  `agglomerative_cluster()`.
- **TRexSelector** (1.0.0, CRAN) — the old package with the functional API
  (`trex()`, `FDP()`, `TPP()`, `add_dummies_GVS()`). The legacy demo areas
  (see table below) are still written against this package; their migration to
  TRexSelectorNeo is planned as a later phase.

---

## Prerequisites

Install **TRexSelectorNeo** from the main repository:

```bash
R CMD INSTALL <path-to-TRexSelector>/R/TRexSelectorNeo
# or, from R:
# devtools::install("<path-to-TRexSelector>/R/TRexSelectorNeo")
```

Additional CRAN packages, depending on which demos you run:

| Package | Needed by |
|---|---|
| `TRexSelector` (1.0.0) | Legacy areas: `trex_selector_methods/trex_da/`, `trex_gvs/`, `trex_spca/`, and `tsolvers/` |
| `plotly` | Correlation-matrix heatmaps in the legacy trex_da / trex_gvs demos |
| `parallel` | Monte Carlo runners (shipped with base R) |
| `glmnet` | `ml_methods/model_selection/` and the `trex_spca/` cross-check probes |
| `tlars` | `tsolvers/demo_ts_compare_tlars_tlasso.R` (independent reference implementation) |
| `elasticnet` | `trex_selector_methods/trex_spca/demo_trex_spca_01.R` (SPCA baseline) |

---

## Running the demos

Every demo is a standalone script and can be run from **any** working directory:

```bash
Rscript R/trex_selector_methods/trex/demo_trex_01_single_run.R
```

Each script resolves its own location (via `sys.frame(1L)$ofile` with a
`--file=` fallback), so relative `source()` calls and output paths work under
both `Rscript` and interactive `source()`. Monte Carlo demos write their
results into a `simulation_results/` folder next to the demo file.

Note on indexing: all supports and selected variables are **1-based** in R;
the C++ and Python counterparts are 0-based. Demo headers state the conversion
explicitly where fixed supports are mirrored from C++.

---

## Contents

| Area | Status | What it covers | cpp counterpart |
|---|---|---|---|
| [memory_mapping/](memory_mapping/README.md) | NEW (TRexSelectorNeo) | `mmap_matrix` / `convert_to_memory_mapped` walkthrough, out-of-core fill, element access | [cpp/memory_mapping/](../cpp/memory_mapping/) |
| [ml_methods/hac_clustering/](ml_methods/hac_clustering/README.md) | NEW (TRexSelectorNeo) | Agglomerative hierarchical clustering (in-memory + mmap), ARI evaluation | [cpp/ml_methods/hac_clustering/](../cpp/ml_methods/hac_clustering/) |
| [ml_methods/standardization/](ml_methods/standardization/README.md) | NEW (TRexSelectorNeo) | `ZScoreScaler` / `LpNormScaler` (in-memory + mmap) | [cpp/ml_methods/scaler_methods/](../cpp/ml_methods/scaler_methods/) |
| [ml_methods/pca/](ml_methods/pca/README.md) | NEW (TRexSelectorNeo) | `PCA` R6 class: fit, restore round-trip, out-of-sample transform | [cpp/ml_methods/pca/](../cpp/ml_methods/pca/) |
| [ml_methods/svd/](ml_methods/svd/README.md) | NEW (TRexSelectorNeo) | Truncated SVD via `SVDSolver` | [cpp/ml_methods/svd/](../cpp/ml_methods/svd/) |
| [ml_methods/ridge_regression/](ml_methods/ridge_regression/README.md) | NEW (TRexSelectorNeo) | L2-penalized least squares via `RidgeSolver` | [cpp/ml_methods/ridge_regression/](../cpp/ml_methods/ridge_regression/) |
| [ml_methods/model_selection/](ml_methods/model_selection/README.md) | NEW (TRexSelectorNeo) + cross-check (glmnet) | K-fold CV ridge via `RidgeCV` (lambda.min / lambda.1se); glmnet reference dumps for the C++ CCD elastic net | [cpp/ml_methods/model_selection/](../cpp/ml_methods/model_selection/) |
| [trex_selector_methods/trex/](trex_selector_methods/trex/README.md) | NEW (TRexSelectorNeo) | Classical T-Rex: single runs, MC studies, l-loop strategies, mmap, scalability (demos 00–07) | [cpp/trex_selector_methods/trex/](../cpp/trex_selector_methods/trex/) |
| [trex_selector_methods/trex_da/](trex_selector_methods/trex_da/README.md) | Legacy (CRAN TRexSelector) | Dependency-Aware T-Rex: AR(1), equi, NN, BT-block, heavy-tailed DGPs (demos 01–08) | [cpp/trex_selector_methods/trex_da/](../cpp/trex_selector_methods/trex_da/) |
| [trex_selector_methods/trex_gvs/](trex_selector_methods/trex_gvs/README.md) | Legacy (CRAN TRexSelector) | Grouped Variable Selection T-Rex: EN vs IEN over grouped DGPs (demos 01–08) | [cpp/trex_selector_methods/trex_gvs/](../cpp/trex_selector_methods/trex_gvs/) |
| [trex_selector_methods/trex_spca/](trex_selector_methods/trex_spca/README.md) | Legacy (CRAN TRexSelector) | T-Rex Sparse PCA MC study + C++/R cross-check probes | [cpp/trex_selector_methods/trex_spca/](../cpp/trex_selector_methods/trex_spca/), [cpp/trex_selector_methods/validation/trex_spca/](../cpp/trex_selector_methods/validation/trex_spca/) |
| [tsolvers/](tsolvers/README.md) | Cross-check (CRAN tlars) | TLARS / TLASSO / TENET path-equivalence reference generator vs C++ tsolvers | [cpp/tsolvers/](../cpp/tsolvers/) |

Shared helpers for the selector demos live in
[trex_selector_methods/](trex_selector_methods/README.md)
(`support_generators.R`, `simulation_utils.R`).

---

## Coverage relative to the C++ tree

The R tree does not yet mirror everything in [cpp/](../cpp/README.md):

- **trex_screening** — the folder `trex_selector_methods/trex_screening/`
  exists but is empty; the six C++ screening demos have no R counterpart yet.
- **tsolvers** — C++ has 12 per-solver demos (`demo_ts_01` … `demo_ts_12`);
  R currently ships only the single TLARS/TLASSO/TENET cross-comparison script.
- **ml_methods** is fully mirrored (hac_clustering, standardization, pca, svd,
  ridge_regression, model_selection). The R elastic-net demo covers K-fold CV
  (`ElasticNetCV`) only; the standalone elastic-net *path* fit from
  `demo_mlm_ms_02_enet_cv_ccd` has no R class and is shown in the Python port
  instead.

---

**Last updated**: 2026-07-06
