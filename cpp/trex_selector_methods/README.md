# T-Rex Selector Methods: Demos and Validation

## Overview

This folder contains all **T-Rex selector variants** built on top of the `tsolvers` and `ml_methods` layers of the C++ project.

Each variant addresses a different variable-selection scenario — the classical baseline, dependent/correlated designs, grouped-variable structures, ultra-high-dimensional screening, and sparse PCA — while sharing the same core idea: control the false discovery rate (FDR) using the Terminating Random Experiments (T-Rex) framework. Demo programs illustrate usage and typical results for each variant, while validation programs cross-check C++ behavior against R reference implementations.

---

## Category overview

| Category | Folder | Purpose |
|----------|--------|---------|
| **Classical T-Rex** | [trex/](trex/README.md) | Terminating Random Experiments Selector — baseline FDR-controlled variable selection in low- and high-dimensional Gaussian linear models |
| **Dependency-Aware T-Rex (DA-TRex)** | [trex_da/](trex_da/README.md) | Variable selection under correlated/dependent design structures: AR(1) Toeplitz, equicorrelated, banded nearest-neighbor, and hierarchical block-covariance designs |
| **Grouped Variable Selection (T-Rex+GVS)** | [trex_gvs/](trex_gvs/README.md) | Group-structured selection via elastic-net-based and HAC-clustering-based grouping (EN, EN-AUG, IEN) across equicorrelated, scattered, and block-covariance patterns |
| **Screening (Screen-TRex)** | [trex_screening/](trex_screening/README.md) | Ultra-high-dimensional pre-screening ahead of T-Rex selection, with in-memory, memory-mapped, and biobank-scale variants |
| **Sparse PCA (T-Rex SPCA)** | [trex_spca/](trex_spca/README.md) | Sparse PCA via T-Rex+GVS applied to loading matrices, benchmarked against ordinary PCA and oracle-thresholded baselines |
| **Validation** | [validation/](validation/README.md) | Cross-checks of C++ selector variants against R reference outputs (currently covers `trex`, `trex_da`, `trex_gvs`, and `trex_spca`) |

---

## Architecture context

These components form the `trex_selector_methods` layer in the project structure:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- each selector variant dispatches to one or more base solvers from `tsolvers` (e.g. TLARS, TLASSO, TENET) within its T-loop,
- variants that rely on preprocessing reuse components from `ml_methods`, for example HAC clustering for group construction in `trex_gvs`, or scaling utilities shared across variants,
- `trex_screening` and the memory-mapped demos in `trex` build on the memory-mapping utilities for large design matrices,
- validation programs check C++ selector output against trusted R reference pipelines.

---

## Building

```bash
cd TRexSelector_Simulations/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

---

## Running

Run any demo or validation executable from the build tree. For example:

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run
```

For each variant's demo lineup and required parameters, see the corresponding subfolder README (once available).

---

## Folder contents

```txt
trex_selector_methods/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex/
  ├── trex_screening/
  ├── trex_da/
  ├── trex_gvs/
  ├── trex_spca/
  └── validation/
```

---

**Last updated**: 2026-07-04
