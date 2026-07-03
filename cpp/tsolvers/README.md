# T-Solvers: Demos and Validation

## Overview

This folder contains demos and validation programs for the **terminating solvers** in `tsolvers/`.

These are variable-selection algorithms designed for high-dimensional linear-model problems inside the T-Rex framework. Each solver follows a different path-following, greedy, or pursuit-style strategy, but all share the same terminating principle: instead of computing a full regularization path, they stop once a prescribed stopping level $T_{\text{stop}}$ has been reached.

The demos in this folder run each solver in a standalone way, outside the full T-Rex calibration workflow, so that solver behavior, normalization choices, and serialization can be inspected in isolation.

---

## Demo overview

| # | Solver | Folder | Strategy |
|---|--------|--------|----------|
| 01 | **T-LARS** | [demo_ts_01_tlars/](demo_ts_01_tlars/README.md) | Least Angle Regression |
| 02 | **T-LASSO** | [demo_ts_02_tlasso/](demo_ts_02_tlasso/README.md) | Lasso path |
| 03 | **T-Stepwise** | [demo_ts_03_tstepwise/](demo_ts_03_tstepwise/README.md) | Forward stepwise regression |
| 04 | **T-ENET** | [demo_ts_04_tenet/](demo_ts_04_tenet/README.md) | Elastic net |
| 05 | **T-Stagewise** | [demo_ts_05_tstagewise/](demo_ts_05_tstagewise/README.md) | Forward stagewise regression |
| 06 | **T-OMP** | [demo_ts_06_tomp/](demo_ts_06_tomp/README.md) | Orthogonal matching pursuit |
| 07 | **T-GP** | [demo_ts_07_tgp/](demo_ts_07_tgp/README.md) | Gradient pursuit |
| 08 | **T-ACGP** | [demo_ts_08_tacgp/](demo_ts_08_tacgp/README.md) | Approximate conjugate gradient pursuit |
| 09 | **T-MP** | [demo_ts_09_tmp/](demo_ts_09_tmp/README.md) | Matching pursuit |
| 10 | **T-NCGMP** | [demo_ts_10_tncgmp/](demo_ts_10_tncgmp/README.md) | Norm-corrective generalized matching pursuit |
| 11 | **T-OLS** | [demo_ts_11_tools/](demo_ts_11_tools/README.md) | Orthogonal least squares |
| 12 | **T-AFS** | [demo_ts_12_tafs/](demo_ts_12_tafs/README.md) | Adaptive forward stepwise |

---

## Validation

In addition to solver demos, this folder also contains heavier validation programs that compare selected C++ solvers against alternative formulations or external reference implementations.

| Folder | Purpose |
|--------|---------|
| [validation/validation_ts_01_tenet_aug_comparison/](validation/validation_ts_01_tenet_aug_comparison/README.md) | Equivalence check between Gram-based `TENET_Solver` and augmented-LASSO `TENETAug_Solver` |
| [validation/validation_ts_02_tlars_tlasso_rcompare/](validation/validation_ts_02_tlars_tlasso_rcompare/README.md) | Cross-language comparison against CRAN `tlars` on sparse-factor-model data |

---

## Key concepts

### Terminating principle

A terminating solver runs its selection path until a stopping level $T_{\text{stop}}$ is reached, which stops it
usually before the full path is computed.
The actual coefficient vector is computed internally, but serves only as tool for the selection process.

### Normalization modes

The demos are organized around two common usage patterns:

- **External normalization**: the caller preprocesses $\mathbf{X}$ before passing it to the solver.
- **Internal normalization**: the solver handles normalization itself and may expose both normalized and raw coefficient views.

### Serialization

Several solver demos also illustrate Cereal-based serialization workflows such as `save(...)` and reload operations, which are relevant for checkpointing and repeated solver use inside larger pipelines.

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
./build/debug/bin/tsolvers/demo_ts_01_tlars/demo_ts_01_tlars
```

For validation executables, see the corresponding README files for any required reference-data inputs.

---

## Folder contents

```txt
tsolvers/
  ├── README.md
  ├── CMakeLists.txt
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
  ├── demo_ts_12_tafs/
  └── validation/
```

---

**Last updated**: 2026-07-03
