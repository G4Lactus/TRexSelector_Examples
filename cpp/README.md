# TRexSelector_Examples — C++ Sources

## Overview

This directory contains the standalone C++ examples project. It consumes the
installed **TRexSelector** library via CMake's `find_package()` and covers the
full library stack: data-preprocessing utilities, T-Algorithm solvers, and T-Rex
selector variants. Each demo writes its output to a `simulation_results/` folder
next to its own source file.

---

## Prerequisites

Install the TRexSelector library before configuring this project:

```bash
cmake --build <TRexSelector-build-dir> --target install
```

The CMake presets in this directory expect the install tree at
`../../TRexSelector/cpp/install`. Adjust `CMAKE_PREFIX_PATH` in
`CMakePresets.json` if your TRexSelector installation is elsewhere.

---

## Build System

### Presets (`CMakePresets.json`)

| Configure preset | Build preset | Build type | Binary output directory |
|---|---|---|---|
| `debug` | `debug-build` | Debug (`-g -O0`) | `build/debug/bin/` |
| `release` | `release-build` | Release (`-O3 -march=native`) | `build/release/bin/` |

Quick start (run from this `cpp/` directory):

```bash
# Full debug build
cmake --preset debug
cmake --build --preset debug-build

# Full release build
cmake --preset release
cmake --build --preset release-build
```

The VS Code tasks **CMake Full Build (Debug/Release)** in
`.vscode/tasks.json` run the same two-step configure-then-build sequence.

### `add_demo_target` macro

All demo executables are registered through the
`add_demo_target(target_name source_file library_set)` macro in
[CMakeLists.txt](CMakeLists.txt). It:

- Creates the executable and links Eigen3, OpenMP, and a chosen library set
- Routes each binary into a mirrored sub-tree under `build/<type>/bin/`
  (e.g. `tsolvers/demo_ts_01_tlars/demo_ts_01_tlars`)
- Injects `DEMO_OUTPUT_DIR` as a compile-time constant pointing to the
  source-side `simulation_results/` folder next to the demo source

### Library sets

| Variable | Linked TRexSelector components |
|---|---|
| `DATA_PREPROCESSING_LIBS` | `trex_ml_methods` |
| `T_ALGORITHM_LIBS` | `trex_ml_methods`, `trex_tsolvers` |
| `T_REX_SELECTOR_LIBS` | `trex_selector_methods` (transitively pulls in tsolvers, ml\_methods, utils) |

---

## Module Structure

### `memory_mapping/`

Memory-mapped matrix demos using disk-backed `Eigen::Map`. Demonstrates
out-of-core matrix creation, zero-copy read/write access, and OpenMP parallel
filling via POSIX `mmap`.

| Demo | Description |
|---|---|
| `demo_memory_mapping` | Create, access, parallel-fill, and clean up a memory-mapped matrix |

---

### `ml_methods/`

Data-preprocessing and numerical building blocks used by the solver and
selector layers.

| Demo | Description |
|---|---|
| `demo_mlm_scaler_01` | Column normalization — Z-score and Lp-norm scalers |
| `demo_mlm_pca_01` | Principal component analysis (PCA) |
| `demo_mlm_svd_01` | Singular value decomposition (SVD) |
| `demo_mlm_ridge_01` | Ridge regression (L2-penalized least squares) |
| `demo_mlm_ms_01_ridge_cv_svd` | Cross-validated ridge regression via SVD |
| `demo_mlm_ms_02_enet_cv_ccd` | Cross-validated elastic net via CCD |
| `demo_mlm_hac_01` | Hierarchical agglomerative clustering (HAC) |
| `demo_mlm_hac_02_mmap` | HAC on a memory-mapped design matrix |

Validation programs in `ml_methods/validation/` cross-check HAC output
against R reference runs.

---

### `tsolvers/`

Terminating T-Algorithm solvers for high-dimensional linear-model variable
selection. Each solver stops at $T_{\text{stop}}$ level rather than computing
the full regularization path.

| Demo | Solver |
|---|---|
| `demo_ts_01_tlars` | T-LARS (Least Angle Regression) |
| `demo_ts_02_tlasso` | T-LASSO |
| `demo_ts_03_tstepwise` | T-Stepwise (forward stepwise) |
| `demo_ts_04_tenet` | T-ENET (elastic net) |
| `demo_ts_05_tstagewise` | T-Stagewise (forward stagewise) |
| `demo_ts_06_tomp` | T-OMP (orthogonal matching pursuit) |
| `demo_ts_07_tgp` | T-GP (gradient pursuit) |
| `demo_ts_08_tacgp` | T-ACGP (approximate conjugate gradient pursuit) |
| `demo_ts_09_tmp` | T-MP (matching pursuit) |
| `demo_ts_10_tncgmp` | T-NCGMP (norm-corrective generalized matching pursuit) |
| `demo_ts_11_tools` | T-OLS (orthogonal least squares) |
| `demo_ts_12_tafs` | T-AFS (adaptive forward stepwise) |

Validation programs in `tsolvers/validation/` check solver equivalences and
cross-compare output against R.

---

### `trex_selector_methods/`

Five T-Rex selector variants, each in its own sub-directory. Demos run Monte
Carlo simulations and write results to source-side `simulation_results/` folders.

#### `trex/`

Classical T-Rex selector. Demos cover single runs, Monte Carlo simulations
with fixed and variable support, l-loop strategies, and memory-mapped designs
(7 demos, 1 pending).

#### `trex_da/`

Dependency-Aware T-Rex. Demos explore AR(1), equicorrelated, nearest-neighbor,
block-diagonal, and group correlation structures (8 demos).

#### `trex_gvs/`

Grouped Variable Selection T-Rex. Demos benchmark EN, EN-AUG, and IEN grouping
strategies across Hastie, scattered, mixed, negative-trap, t-distributed, AR(1),
ARMA, and block-benchmark designs (8 demos).

#### `trex_screening/`

Ultra-high-dimensional screening with T-Rex. Demos cover in-memory screening,
memory-mapped screening, correlated designs, and biobank-scale pipelines
(6 demos).

#### `trex_spca/`

Sparse PCA with T-Rex. Monte Carlo simulation demo (1 demo).

Validation programs in `trex_selector_methods/validation/` cover scaling
comparisons, solver equivalence checks, and R cross-validation pipelines.

---

*Last updated: 2026-07-04*
