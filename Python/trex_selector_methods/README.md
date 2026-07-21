# T-Rex Selector Methods: Demos

## Overview

This folder contains the Python demos for all **T-Rex selector variants**, built on the `trex_selector_neo` package (the pybind11 binding layer over the C++ `tsolvers` and `ml_methods` stack).

Each variant addresses a different variable-selection scenario — the classical baseline, dependent/correlated designs, grouped-variable structures, ultra-high-dimensional screening, and sparse PCA — while sharing the same core idea: control the false discovery rate (FDR) using the Terminating Random Experiments (T-Rex) framework. The demos mirror their C++ counterparts in [cpp/trex_selector_methods/](../../cpp/trex_selector_methods/): same data-generating processes, same control parameters, same result summaries, following the same **folder-per-demo** layout (one subfolder per demo, each with its own `README.md` and `simulation_results/`, plus suite-level shared modules).

The C++ tree's `validation/` programs (cross-checks against R reference outputs) are validation infrastructure rather than demos and have no Python port by design; they live in the TRexSelector library test suite.

---

## Category overview

| Category | Folder | Purpose |
|----------|--------|---------|
| **Classical T-Rex** | [trex/](trex/README.md) | Terminating Random Experiments Selector — baseline FDR-controlled variable selection in low- and high-dimensional Gaussian linear models (demos 01–08) |
| **Dependency-Aware T-Rex (DA-TRex)** | [trex_da/](trex_da/README.md) | Variable selection under correlated/dependent design structures: AR(1) Toeplitz, equicorrelated, banded nearest-neighbor, and hierarchical block-covariance designs (demos 01–08) |
| **Grouped Variable Selection (T-Rex+GVS)** | [trex_gvs/](trex_gvs/README.md) | Group-structured selection via elastic-net-based and HAC-clustering-based grouping (EN, EN-AUG, IEN) across equicorrelated, scattered, and block-covariance patterns (demos 01–08) |
| **Screening (Screen-TRex)** | [trex_screening/](trex_screening/README.md) | Ultra-high-dimensional pre-screening ahead of T-Rex selection, with in-memory, memory-mapped, and biobank-scale variants (demos 01–06) |
| **Sparse PCA (T-Rex SPCA)** | [trex_spca/](trex_spca/README.md) | Sparse PCA via per-component T-Rex+GVS selection on plug-in PCs, benchmarked against ordinary PCA and oracle-thresholded baselines: Fig.-2/Fig.-3-style support-recovery and Definition-1 PEV studies, plus mechanism studies of the plug-in construction (plug-in vs. oracle response, union-support FDR per PC) (demos 01–02) |

---

## Architecture context

The demos exercise the `trex_selector_neo` selector layer, which sits on top of the solver and preprocessing layers:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- each selector variant dispatches to one or more base solvers from `tsolvers` (e.g. TLARS, TLASSO, TENET) within its T-loop,
- variants that rely on preprocessing reuse components from `ml_methods`, for example HAC clustering for group construction in `trex_gvs`, or scaling utilities shared across variants,
- `trex_screening` and the memory-mapped demos in `trex` build on the memory-mapping utilities for large design matrices.

All selector-layer classes and enums are accessed through the package root (`import trex_selector_neo as tsn`); see the [import conventions](../README.md#import-conventions).

---

## Running

Each demo is a plain script, run directly with the Python interpreter:

```bash
python trex/demo_trex_01_single_run/demo_trex_01_single_run.py
```

For each variant's demo lineup and required parameters, see the corresponding subfolder README.

---

## Folder contents

```txt
trex_selector_methods/
  ├── README.md
  ├── trex/
  ├── trex_screening/
  ├── trex_da/
  ├── trex_gvs/
  └── trex_spca/
```

---

**Last updated**: 2026-07-21
