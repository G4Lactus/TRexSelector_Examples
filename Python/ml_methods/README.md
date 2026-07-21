# ML Methods: Demos

## Overview

This folder contains the Python demos for the `ml_methods` layer of the TRexSelector library, exercised through the
`trex_selector_neo.ml_methods` module.

It provides the numerical and preprocessing components that support higher-level T-Rex workflows, including column
scaling, singular value decomposition, principal component analysis, ridge regression, model parameter selection based
on k-fold cross-validation, and hierarchical agglomerative clustering.

Each method family has its own subfolder. Demo programs illustrate usage and interface behavior. Cross-language
correctness checks against R references live in the TRexSelector library's own test suite (`cpp/tests/validation/`),
run on demand via its `run_validations` target, rather than in this examples package.

The demos mirror the C++ suite in [cpp/ml_methods/](../../cpp/ml_methods/): same sub-demo structure, parameters, and
console-output layout, with Eigen matrices replaced by Fortran-ordered NumPy arrays. Where the C++ demos run at very
large scale (e.g., OpenMP-parallel generation of matrices with $p \ge 50{,}000$), the Python mirrors scale the problem
sizes down — the reductions are noted in each script's header comment.

---

## Category overview

| Category | Folder | Purpose |
| -------- | ------ | ------- |
| **Scaler Methods** | [scaler_methods/](scaler_methods/README.md) | Column normalization utilities, including Z-score and Lp-norm scaling |
| **PCA** | [pca/](pca/README.md) | Principal component analysis |
| **SVD** | [svd/](svd/README.md) | Singular value decomposition for stable matrix computations |
| **Ridge Regression** | [ridge_regression/](ridge_regression/README.md) | L2-penalized least-squares regression |
| **Model Selection** | [model_selection/](model_selection/README.md) | Cross-validation workflows for regularized regression methods |
| **HAC Clustering** | [hac_clustering/](hac_clustering/README.md) | Hierarchical agglomerative clustering, including memory-mapped workflows |

---

## Architecture context

These components form the `ml_methods` layer in the project structure:

```txt
utils -> ml_methods -> tsolvers -> trex_selector_methods -> demos
```

In particular:

- scaler methods support feature normalization before downstream modeling,
- SVD and ridge-related tools support stable regression and dimensionality-reduction workflows,
- HAC clustering supports grouped-variable workflows such as GVS.

---

## Running

Each demo is a plain script run directly with the Python interpreter (no build step required), e.g.:

```bash
python pca/demo_mlm_pca_01/demo_mlm_pca_01.py
```

All demos assume an environment with `trex_selector_neo` (v0.2.0) and `numpy` installed; see the
[Python README](../README.md) for installation and the mandatory import conventions.

---

## Folder contents

```txt
ml_methods/
  ├── README.md
  ├── scaler_methods/
  ├── pca/
  ├── svd/
  ├── ridge_regression/
  ├── model_selection/
  └── hac_clustering/
```

---

**Last updated**: 2026-07-21
