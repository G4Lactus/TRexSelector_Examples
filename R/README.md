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
  (`trex()`, `FDP()`, `TPP()`, `add_dummies_GVS()`). Its only consumer, the
  path-equivalence reference generator `demo_ts_compare_tlars_tlasso.R`, moved to
  the TRexSelector library test suite
  (`TRexSelector/cpp/tests/validation/tsolvers/`), so no demo in this examples
  tree depends on it any more. All the selector suites (trex_da, trex_gvs,
  trex_spca, and the new trex_screening) have been migrated to the
  TRexSelectorNeo R6 API and no longer depend on it.

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
| `TRexSelector` (1.0.0) | the reference generator `demo_ts_compare_tlars_tlasso.R`, now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`) |
| `plotly` | Correlation-matrix heatmaps in the trex_da / trex_gvs MC demos |
| `parallel` | Monte Carlo runners (shipped with base R) |
| `glmnet` | the model-selection and trex_spca cross-check probes, now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/ml_methods/model_selection/` and `.../trex_selector_methods/trex_spca/`) |
| `tlars` | the reference generator `demo_ts_compare_tlars_tlasso.R`, now in the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`) |

---

## Running the demos

Every demo is a standalone script and can be run from **any** working directory:

```bash
Rscript R/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run.R
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
| [trex_selector_methods/trex_da/](trex_selector_methods/trex_da/README.md) | NEW (TRexSelectorNeo) | Dependency-Aware T-Rex: AR(1), equi, NN, BT-block, heavy-tailed DGPs (demos 01–08); NN demos skip until the R binding exposes `da_method="NN"` | [cpp/trex_selector_methods/trex_da/](../cpp/trex_selector_methods/trex_da/) |
| [trex_selector_methods/trex_gvs/](trex_selector_methods/trex_gvs/README.md) | NEW (TRexSelectorNeo) | Grouped Variable Selection T-Rex: EN vs IEN over grouped DGPs (demos 01–08) | [cpp/trex_selector_methods/trex_gvs/](../cpp/trex_selector_methods/trex_gvs/) |
| [trex_selector_methods/trex_spca/](trex_selector_methods/trex_spca/README.md) | NEW (TRexSelectorNeo) | T-Rex Sparse PCA MC comparison vs. ordinary/oracle PCA (`TRexSPCASelector`, demo 01, folder-per-demo); C++/R cross-check probes moved to the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`) | [cpp/trex_selector_methods/trex_spca/](../cpp/trex_selector_methods/trex_spca/) |
| [trex_selector_methods/trex_screening/](trex_selector_methods/trex_screening/README.md) | NEW (TRexSelectorNeo) | Screen-TRex: in-memory / mmap screening, correlated designs, biobank-scale, solver comparison (6 demos) | [cpp/trex_selector_methods/trex_screening/](../cpp/trex_selector_methods/trex_screening/) |
| [tsolvers/](tsolvers/README.md) | NEW (TRexSelectorNeo) + cross-check (CRAN tlars) | 12 standalone terminating-solver demos (`demo_ts_01` … `demo_ts_12`: early stopping, external normalization, serialization, mmap); the TENET/TENETAug equivalence check and CRAN tlars path-equivalence reference generator moved to the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`) | [cpp/tsolvers/](../cpp/tsolvers/) |

Shared helpers for the selector demos live in
[trex_selector_methods/](trex_selector_methods/README.md)
(`support_generators.R`, `simulation_utils.R`).

---

## Coverage relative to the C++ tree

The R tree now mirrors the selector suites — `trex`, `trex_da`, `trex_gvs`,
`trex_spca`, and `trex_screening` — on the new TRexSelectorNeo R6 API. A few
gaps relative to [cpp/](../cpp/README.md) remain:

- **trex+DA+NN** — the C++ core supports `DAMethod::NN`, but the installed
  TRexSelectorNeo R binding does not expose `da_method="NN"` yet, so the two
  `trex_da` nearest-neighbour demos skip cleanly until it is added.
- **tsolvers** — all 12 per-solver demos are mirrored
  (`demo_ts_01` … `demo_ts_12`). The validation programs
  (`validation_ts_01_tenet_aug_comparison` and the CRAN tlars cross-comparison
  generator `validation_ts_02_tlars_tlasso_rcompare`) moved to the TRexSelector
  library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`).
- **ml_methods** is fully mirrored (hac_clustering, standardization, pca, svd,
  ridge_regression, model_selection); pca / svd / ridge are covered there. The
  R elastic-net demo covers K-fold CV (`ElasticNetCV`) only; the standalone
  elastic-net *path* fit from `demo_mlm_ms_02_enet_cv_ccd` has no R class and is
  shown in the Python port instead.

---

**Last updated**: 2026-07-08
