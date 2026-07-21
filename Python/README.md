# TRexSelector_Examples — Python Sources

## Overview

This directory contains the standalone Python examples project. It consumes the
**TRexSelectorNeo** package (v0.2.0), the pybind11 binding layer over the C++
TRexSelector backend, imported as

```python
import trex_selector_neo
```

The demos mirror their C++ counterparts in the sibling [cpp/](../cpp/) tree
(and the R demos in [R/](../R/)): same demo names and folder-per-demo layout
(one subfolder per demo, each with its own `README.md`), same data-generating
processes, same control parameters, same result summaries. All demos are plain
scripts — no build step is required. Monte Carlo demos write their output to a
`simulation_results/` folder next to the demo script.

> **Indexing note**: All variable indices in the Python demos are **0-based**
> (matching C++). The R package is 1-based, so R support sets are shifted by
> one relative to the values printed here.

---

## Import Conventions

All demos follow one canonical import scheme so that every symbol's origin is
visible at its use site:

```python
import os
import sys

import numpy as np
import pandas as pd

import trex_selector_neo as tsn
from trex_selector_neo.ml_methods import LpNormScaler, NormType
from trex_selector_neo.ml_methods.clustering import LinkageMethod, agglomerative_cluster
from trex_selector_neo.tsolvers.lars_based import TLARS_Solver
from trex_selector_neo.utils import compute_fdp, compute_tpp, numpy_to_memmap

from trex_sim_common import run_mc_trex  # local demo helpers last
```

The rules:

1. **Selector layer via the `tsn` alias.** Everything the package root exports —
   selector classes (`tsn.TRexSelector`, `tsn.TRexDASelector`, …), control
   parameters (`tsn.TRexControlParameter`, …), and enums
   (`tsn.SolverTypeForTRex`, `tsn.LLoopStrategy`, `tsn.DummyDistribution`, …) —
   is accessed through `import trex_selector_neo as tsn`, mirroring the C++
   `namespace` qualification and R's `trexselector::`. Root symbols are never
   `from`-imported, and the alias is only introduced in files that use it.
2. **Subpackage symbols via explicit `from`-imports.** Solvers come from
   `trex_selector_neo.tsolvers.{lars_based,omp_based,afs_based}`, ML building
   blocks from `trex_selector_neo.ml_methods` (clustering from
   `trex_selector_neo.ml_methods.clustering`), and metrics/memory-mapping
   helpers from `trex_selector_neo.utils`. Call sites use the bare names;
   subpackages are never aliased and never accessed as `tsn.utils.…`.
3. **Only the public surface.** Demos never import from the raw binding module
   `trex_selector_neo.trex_selector_methods` (it bypasses the ergonomic Python
   wrappers), never from `trex_selector_neo._core`, and never through the
   legacy `ml_methods.clustering.hierarchical` shim.
4. **Import order.** Standard library → third-party (`numpy`, `pandas`,
   `matplotlib`) → `trex_selector_neo` (root alias first, then subpackages
   alphabetically) → local demo helper modules.

---

## Prerequisites and Installation

- Python >= 3.10
- numpy >= 1.26 (installed automatically with TRexSelectorNeo)
- pandas (used by the demos to write CSV result files)

TRexSelectorNeo is built with scikit-build-core + pybind11 from the main
repository. Install it into a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install <path-to>/TRexSelector/Python/TRexSelectorNeo
pip install pandas
```

The old package name `trex_selector` is obsolete — all demos in this tree
import `trex_selector_neo`.

---

## Running a Demo

Each demo is run directly with the Python interpreter, e.g.:

```bash
python trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run.py
```

Notes:

- Each demo inserts its own directory into `sys.path`, so sibling helper
  modules resolve regardless of the current working directory.
- The Monte Carlo simulation demos parallelize trials with
  `concurrent.futures.ProcessPoolExecutor` and write results
  (pandas CSV + plain-text tables) to a `simulation_results/` folder
  next to the demo script.

---

## Contents

| Area | Folder | What it covers | C++ counterpart |
|---|---|---|---|
| Memory mapping | [memory_mapping/](memory_mapping/) | Disk-backed matrices via `MemoryMappedMatrix`, `numpy_to_memmap`, out-of-core workflows | [cpp/memory_mapping/](../cpp/memory_mapping/) |
| ML methods — HAC | [ml_methods/hac_clustering/](ml_methods/hac_clustering/) | Hierarchical agglomerative clustering: sample/variable clustering, LSH linkage, linkage comparisons, out-of-core mmap workflow (`demo_mlm_hac_01`, `demo_mlm_hac_02_mmap`) | [cpp/ml_methods/hac_clustering/](../cpp/ml_methods/hac_clustering/) |
| ML methods — Scalers | [ml_methods/scaler_methods/](ml_methods/scaler_methods/) | Z-score and Lp-norm column scaling, inverse transform, serialization, mmap-backed (`demo_mlm_scaler_01`) | [cpp/ml_methods/scaler_methods/](../cpp/ml_methods/scaler_methods/) |
| ML methods — PCA | [ml_methods/pca/](ml_methods/pca/) | PCA fit, restore round-trip, out-of-sample transform (`demo_mlm_pca_01`) | [cpp/ml_methods/pca/](../cpp/ml_methods/pca/) |
| ML methods — SVD | [ml_methods/svd/](ml_methods/svd/) | Truncated SVD, direct and Gram paths (`demo_mlm_svd_01`) | [cpp/ml_methods/svd/](../cpp/ml_methods/svd/) |
| ML methods — Ridge | [ml_methods/ridge_regression/](ml_methods/ridge_regression/) | L2-penalized least squares, primal/dual, lambda sweep (`demo_mlm_ridge_01`) | [cpp/ml_methods/ridge_regression/](../cpp/ml_methods/ridge_regression/) |
| ML methods — Model selection | [ml_methods/model_selection/](ml_methods/model_selection/) | K-fold CV for ridge (SVD) and elastic net (CCD), lambda.min / lambda.1se (`demo_mlm_ms_01`, `demo_mlm_ms_02`) | [cpp/ml_methods/model_selection/](../cpp/ml_methods/model_selection/) |
| T-Rex selector | [trex_selector_methods/trex/](trex_selector_methods/trex/) | Classical T-Rex: single runs, Monte Carlo simulations, L-loop strategies, dummy distributions, memory-mapped pipelines, scalability benchmark (Demos 01–08) | [cpp/trex_selector_methods/trex/](../cpp/trex_selector_methods/trex/) |
| Dependency-Aware T-Rex | [trex_selector_methods/trex_da/](trex_selector_methods/trex_da/) | T-Rex+DA: AR(1), AR(1)/heavy-tailed blocks (± white noise), prior groups, nearest-neighbor, equicorrelated + BT (Demos 01–08) | [cpp/trex_selector_methods/trex_da/](../cpp/trex_selector_methods/trex_da/) |
| Grouped Variable Selection T-Rex | [trex_selector_methods/trex_gvs/](trex_selector_methods/trex_gvs/) | T-Rex+GVS: Hastie, scattered, mixed/negative-trap, heavy-tailed, AR(1), ARMA, block benchmark; EN / EN+AUG / IEN (Demos 01–08) | [cpp/trex_selector_methods/trex_gvs/](../cpp/trex_selector_methods/trex_gvs/) |
| Screen-TRex | [trex_selector_methods/trex_screening/](trex_selector_methods/trex_screening/) | Ultra-high-dimensional screening: Ordinary/Bootstrap-CI, mmap, correlated designs, biobank routing, solver backends (Demos 01–06) | [cpp/trex_selector_methods/trex_screening/](../cpp/trex_selector_methods/trex_screening/) |
| T-Rex Sparse PCA | [trex_selector_methods/trex_spca/](trex_selector_methods/trex_spca/) | T-Rex SPCA: MC comparison vs. ordinary/oracle PCA on a sparse 3-factor model, plus Fig.-3-style PEV sweeps and the plugin-vs-oracle FDR/TPR heatmap study (Demos 01–02) | [cpp/trex_selector_methods/trex_spca/](../cpp/trex_selector_methods/trex_spca/) |
| T-Solvers | [tsolvers/](tsolvers/README.md) | 12 standalone terminating-solver demos (`demo_ts_01` … `demo_ts_12`: early stopping, external normalization, serialization, mmap); the validation programs moved to the TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`) | [cpp/tsolvers/](../cpp/tsolvers/) |

---

## Coverage vs. cpp

The Python examples cover the classical T-Rex selector, the Dependency-Aware
T-Rex, the Grouped-Variable-Selection T-Rex, Screen-TRex, T-Rex Sparse PCA,
the standalone terminating solvers (`tsolvers/`, 12 demos), memory mapping,
and the full `ml_methods` demo set (HAC clustering, scaler methods, PCA, SVD,
ridge regression, model selection) — every cpp demo suite is ported, with
matching demo names, folder layout, and result-file schemas (as of the
2026-07-21 sync).

The scaler/HAC cross-language reference programs that used to live under
`cpp/ml_methods/validation/` moved to the TRexSelector library test suite
(`TRexSelector/cpp/tests/validation/hac_clustering/`); as validation
infrastructure rather than demos, they have no Python port by design.

---

**Last updated**: 2026-07-21
