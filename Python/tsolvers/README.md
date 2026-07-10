# T-Solvers: Demos (Python)

## Overview

This folder contains the Python demos for the
**terminating solvers** exposed by `trex_selector_neo.tsolvers`. It mirrors
[cpp/tsolvers/](../../cpp/tsolvers/README.md) 1:1 (folder-per-demo, same
numbering and content).

These are variable-selection algorithms designed for high-dimensional
linear-model problems inside the T-Rex framework. Each solver follows a
different path-following, greedy, or pursuit-style strategy, but all share the
same terminating principle: instead of computing a full regularization path,
they stop early once a prescribed number $T_{\text{stop}}$ of dummy variables
has entered the active set.

The demos run each solver in a standalone way, outside the full T-Rex
calibration workflow, so that solver behavior, normalization choices, and
serialization can be inspected in isolation.

---

## Demo overview

| # | Solver | Folder | Strategy | Binding class |
|---|--------|--------|----------|---------------|
| 01 | **T-LARS** | [demo_ts_01_tlars/](demo_ts_01_tlars/README.md) | Least Angle Regression | `tsolvers.lars_based.TLARS_Solver` |
| 02 | **T-LASSO** | [demo_ts_02_tlasso/](demo_ts_02_tlasso/README.md) | Lasso path | `tsolvers.lars_based.TLASSO_Solver` |
| 03 | **T-Stepwise** | [demo_ts_03_tstepwise/](demo_ts_03_tstepwise/README.md) | Forward stepwise regression | `tsolvers.lars_based.TSTEPWISE_Solver` |
| 04 | **T-ENET** | [demo_ts_04_tenet/](demo_ts_04_tenet/README.md) | Elastic net | `tsolvers.lars_based.TENET_Solver` |
| 05 | **T-Stagewise** | [demo_ts_05_tstagewise/](demo_ts_05_tstagewise/README.md) | Forward stagewise regression | `tsolvers.lars_based.TSTAGEWISE_Solver` |
| 06 | **T-OMP** | [demo_ts_06_tomp/](demo_ts_06_tomp/README.md) | Orthogonal matching pursuit | `tsolvers.omp_based.TOMP_Solver` |
| 07 | **T-GP** | [demo_ts_07_tgp/](demo_ts_07_tgp/README.md) | Gradient pursuit | `tsolvers.omp_based.TGP_Solver` |
| 08 | **T-ACGP** | [demo_ts_08_tacgp/](demo_ts_08_tacgp/README.md) | Approximate conjugate gradient pursuit | `tsolvers.omp_based.TACGP_Solver` |
| 09 | **T-MP** | [demo_ts_09_tmp/](demo_ts_09_tmp/README.md) | Matching pursuit | `tsolvers.omp_based.TMP_Solver` |
| 10 | **T-NCGMP** | [demo_ts_10_tncgmp/](demo_ts_10_tncgmp/README.md) | Norm-corrective generalized matching pursuit | `tsolvers.omp_based.TNCGMP_Solver` |
| 11 | **T-OLS** | [demo_ts_11_tools/](demo_ts_11_tools/README.md) | Orthogonal least squares | `tsolvers.omp_based.TOOLS_Solver` |
| 12 | **T-AFS** | [demo_ts_12_tafs/](demo_ts_12_tafs/README.md) | Adaptive forward stepwise | `tsolvers.afs_based.TAFS_Solver` |

Every demo has the same four parts (mirroring the cpp files):

1. **Early stopping** — internal normalization, high-dimensional $(p > n)$ and
   low-dimensional $(n > p)$ regimes.
2. **External normalization** — `LpNormScaler` (center + unit-L2 columns) on
   X and D plus manual y-centering, solver run with
   `normalize=False, intercept=False`.
3. **Serialization & warm start** — save a partial path, reload, continue,
   and compare against a reference run (steps, RSS, $R^2$, action path).
4. **Memory-mapped data** — X/D/y written to disk via `numpy_to_memmap()` and
   solved on the zero-copy mmap views. (cpp generates $p = 500{,}000$ on
   disk; the port scales this to $p = 5{,}000$.)

---

## Validation

The Python validation programs for the terminating solvers moved to the
TRexSelector library test suite (`TRexSelector/cpp/tests/validation/tsolvers/`),
co-located with their R and C++ counterparts:

- `validation_ts_01_tenet_aug_comparison` — equivalence check between Gram-based
  `TENET_Solver`, augmented-LASSO `TENETAug_Solver`, and a manual lasso_star
  augmentation solved with `TLASSO_Solver`.
- `validation_ts_02_tlars_tlasso_rcompare` — path-equivalence gate of the
  Python-binding solvers against CRAN `tlars` reference dumps on
  sparse-factor-model data. The frozen reference dumps now live at
  `TRexSelector/cpp/tests/validation/data/rdump_tlars/`.

---

## Key concepts

### Terminating principle

A terminating solver runs its selection path with
`executeStep(T_stop, early_stop=True)`, where $T_{\text{stop}}$ is the
**number of dummy variables that must enter the active set before the solver
terminates early**, usually well before the full path is computed. Passing
`early_stop=False` traces the complete path instead.

### Shared helpers

[`ts_demo_utils.py`](ts_demo_utils.py) replicates the cpp demo utilities:
`gen_synthetic_data()` (numpy version of `datagen::SyntheticData` with
SNR-calibrated Gaussian noise) and the `print_*` diagnostics
(configuration, selection, and TPP/FDP/precision/recall/F1 quality measures
via `trex_selector_neo.utils`).

### Serialization

The binding serializes solver state with `solver.save(path)` and hydrates an
existing instance with `solver.load(path)` (cpp uses a static
`load(filename, X, D)` factory); the constructor re-binds the data views.

---

## Running

Run any demo from the repository root, e.g.:

```bash
.venv/bin/python Python/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars.py
```

All demos print their diagnostics to the console and clean up any checkpoint
or mmap files they create inside their own folder.

---

## Folder contents

```txt
tsolvers/
  ├── README.md
  ├── ts_demo_utils.py
  ├── demo_ts_01_tlars/
  ├── demo_ts_02_tlasso/
  ├── demo_ts_03_tstepwise/
  ├── demo_ts_04_tenet/
  ├── demo_ts_05_tstagewise/
  ├── demo_ts_06_tomp/
  ├── demo_ts_07_tgp/
  ├── demo_ts_08_tacgp/
  ├── demo_ts_09_tmp/
  ├── demo_ts_10_tncgmp/
  ├── demo_ts_11_tools/
  └── demo_ts_12_tafs/
```

---

**Last updated**: 2026-07-08
