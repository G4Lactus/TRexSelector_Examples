# Model Selection: Demonstration Suite

## Overview

This folder contains C++ examples for **model selection** routines in the `ml_methods` module.

The demos focus on regularized Gaussian regression models and show how hyperparameters are selected from synthetic data
using K-fold cross-validation. In the TRexSelector project, these workflows are relevant for tuning penalization
strength in ridge- and elastic-net-type estimators.

The main goals of this folder are:

1. to demonstrate model selection for regularized linear regression in C++,
2. to illustrate how K-fold cross-validation is used to choose penalty parameters,
3. to compare numerically stable SVD-based ridge fitting with coordinate-descent elastic-net fitting,
4. to connect fitted penalty choices to downstream TRexSelector workflows.

Equivalent model selection demos are also available in the Python and R packages; see
[Python/ml_methods/model_selection/](../../../Python/ml_methods/model_selection/README.md) and
[R/ml_methods/model_selection/](../../../R/ml_methods/model_selection/README.md) for the corresponding walkthroughs.

---

## What this folder covers

The demos study penalized regression models of the form

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n).
$$

They focus on:

- **ridge regression**, where coefficient shrinkage is controlled by a quadratic penalty,
- **elastic-net regression**, which combines $\ell_1$ and $\ell_2$ regularization,
- **cross-validation workflows**, used to select tuning parameters such as $\lambda$ and $\alpha$.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_ms_01_ridge_cv_svd/`**  
   A ridge-regression CV demo using an SVD-based solver across several dimensional regimes.

2. **`demo_mlm_ms_02_enet_cv_ccd/`**  
   An elastic-net demo covering both regularization-path fitting and K-fold CV via coordinate descent.

---

## Demo overview

| Folder | Purpose |
| -------------------------------- | --------- |
| `demo_mlm_ms_01_ridge_cv_svd/` | Ridge regression K-fold CV via `Eigen::JacobiSVD`, including a memory-mapped design-matrix example |
| `demo_mlm_ms_02_enet_cv_ccd/` | Elastic-net path fitting and K-fold CV via coordinate descent, including ridge/lasso/elastic-net comparisons |

---

## Running

```bash
./build/debug/bin/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/demo_mlm_ms_01_ridge_cv_svd
./build/debug/bin/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/demo_mlm_ms_02_enet_cv_ccd
```

---

## Folder contents

```txt
model_selection/
  ├── README.md
  ├── demo_mlm_ms_01_ridge_cv_svd/
  │   ├── demo_mlm_ms_01_ridge_cv_svd.cpp
  │   └── README.md
  └── demo_mlm_ms_02_enet_cv_ccd/
      ├── demo_mlm_ms_02_enet_cv_ccd.cpp
      └── README.md
```

---

## Notes for new users

- Both demos use synthetic Gaussian regression data with known signal structure.
- The focus is on hyperparameter selection and diagnostic output rather than downstream prediction benchmarking.
- The demo-specific README files explain the exact model, penalty, and cross-validation setup in more detail.

---

**Last updated**: 2026-07-01
