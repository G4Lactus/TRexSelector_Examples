# T-Rex Selector Demos and Simulations

This repository contains runnable examples, demos, and Monte Carlo simulation studies for the
[TRexSelector](../TRexSelector) software — a high-performance C++ implementation of the T-Rex
selector family with R and Python bindings. The examples are organized per language and are
intended as the first stop for new users: each demo is a small, self-contained script or program
that generates data, runs a selector (or building block), and explains the output.

---

## Contents

| Folder | Language | Status | Description |
|---|---|---|---|
| [cpp/](cpp/) | C++20 | reference | The most complete example set: all selector variants, 12 terminating solvers, ML building blocks, memory mapping, and validation programs. Has its own CMake build; see [cpp/README.md](cpp/README.md). |
| [R/](R/) | R | partially updated | Demos for the `TRexSelectorNeo` R package. The classical T-Rex, memory-mapping, and ML-method demos use the new R6 API; the `trex_da`/`trex_gvs`/`trex_spca` suites are legacy (CRAN `TRexSelector` 1.0.0) pending migration. See [R/README.md](R/README.md). |
| [Python/](Python/) | Python ≥ 3.10 | updated subset | Demos for the `TRexSelectorNeo` Python package (`import trex_selector_neo`): classical T-Rex, memory mapping, and the full ml_methods set (HAC clustering, standardization, PCA, SVD, ridge regression). Selector-variant areas are not yet ported. See [Python/README.md](Python/README.md). |

Demo naming follows the C++ convention `demo_<area>_NN_<topic>`; R and Python demos state their
C++ counterpart in the file header so the three implementations can be compared side by side.

---

## About the T-Rex Selector

The T-Rex (Terminating-Random EXperiments) selector performs fast, FDR-controlled variable
selection in high-dimensional regression problems (Machkour, Muma, and Palomar, 2022+). It runs
K random experiments in which the design matrix is augmented with computer-generated dummy
variables and a terminating solver (e.g. T-LARS) stops early once enough dummies enter the path;
a calibration procedure then selects a voting threshold that provably controls the false
discovery rate at a user-chosen target level while maximizing the number of selected variables.

The software implements the classical T-Rex selector plus several extensions:

- **T-Rex DA** — dependency-aware variants for correlated designs (AR(1), equicorrelation,
  nearest-neighbor, and binary-tree/dendrogram corrections).
- **T-Rex GVS** — grouped variable selection with elastic-net-type solvers and
  clustering-informed dummies.
- **T-Rex Screening / Biobank Screening** — ultra-high-dimensional screening workflows.
- **T-Rex SPCA** — sparse principal component analysis with FDR control on the loadings.
- **Terminating solvers (tsolvers)** — T-LARS, T-LASSO, T-ENET, matching-pursuit and stepwise
  families, all with early termination at a dummy-count threshold.
- **ML methods** — the numerical building blocks (scalers, PCA, SVD, ridge, elastic-net CV,
  hierarchical agglomerative clustering) and disk-backed memory-mapped matrices for
  larger-than-RAM problems.

---

## Overview of supported Languages and Packages

All three languages are backed by the same C++ core; results agree across languages up to the
documented RNG and indexing conventions (C++/Python are 0-based, R is 1-based).

### R Package

`TRexSelectorNeo` (R ≥ 4.3, Rcpp/RcppEigen, R6 classes). Install from the main repo:

```sh
R CMD INSTALL <path-to>/TRexSelector/R/TRexSelectorNeo
```

Typical usage: `TRexSelector$new(X, y, tFDR = 0.1, control = trex_control(...))$select()`.
Some legacy demos additionally use the CRAN package `TRexSelector` (1.0.0, functional API) —
see [R/README.md](R/README.md) for per-area status.

### Python Package

`TRexSelectorNeo` (Python ≥ 3.10, pybind11/scikit-build-core), imported as `trex_selector_neo`:

```sh
pip install <path-to>/TRexSelector/Python/TRexSelectorNeo
```

Typical usage: `trex_selector_neo.TRexSelector(X, y, tFDR=0.1, trex_control=ctrl).select()`.

### Cpp Package

The C++ library `TRexSelector` is consumed via CMake (`find_package(TRexSelector REQUIRED)`).
The [cpp/](cpp/) folder is a standalone CMake project; build it with the provided presets
(see [cpp/README.md](cpp/README.md)) and run the demo executables from `bin/`.

---

**Last updated**: 2026-07-06
