# Model Selection — R Demos

This folder contains three model-selection items: a **RidgeCV demo** (K-fold
cross-validation for the ridge penalty via the TRexSelectorNeo `RidgeCV` R6
class), an **ElasticNetCV demo** (K-fold CV for the elastic net via
coordinate descent), and an **elastic-net glmnet cross-check** reference
generator.

---

## Ridge K-Fold Cross-Validation (RidgeCV)

Demonstrates ridge-regression K-fold cross-validation with the
**TRexSelectorNeo** `RidgeCV` class (SVD-based solver with per-fold centering
and column L2-normalization): lambda-grid construction, the CV-MSE curve, and
the `lambda.min` / `lambda.1se` selections across three dimensional regimes.

| File | Description |
|---|---|
| [demo_ridge_cv_01.R](demo_ridge_cv_01.R) | K-fold ridge CV on synthetic sparse-signal data — low-dimensional (n=1000, p=500, 10-fold), high-dimensional (n=300, p=500, 10-fold), and very high-dimensional (n=300, p=1000, 5-fold); each scenario prints the CV curve and checks `lambda.min <= lambda.1se` via `stopifnot()` |

### APIs used

- `RidgeCV$new()` with
  `$fit(X, y, num_folds, n_lambda, lambda_ratio, seed)`, `$cv_min()`,
  `$cv_1se()`, `$index_min()`, `$index_1se()`, `$get_lambdas()`,
  `$get_cv_errors()`, `$get_cv_std()`.

The C++ counterpart (`demo_mlm_ms_01_ridge_cv_svd`) additionally runs a
memory-mapped scenario; `RidgeCV`'s R binding only accepts in-memory
matrices, so that scenario is omitted here.

### Running the RidgeCV demo

```bash
Rscript R/ml_methods/model_selection/demo_ridge_cv_01.R
```

Console output only.

---

## Elastic-Net K-Fold Cross-Validation (ElasticNetCV)

Demonstrates elastic-net K-fold cross-validation with the **TRexSelectorNeo**
`ElasticNetCV` class (glmnet-style pathwise cyclic coordinate descent):
CV lambda selection across the alpha range (ridge / elastic net / lasso) and
the pure-ridge `lambda_2_lars = lambda.1se * p / 2` conversion used
internally by TRexGVS.

| File | Description |
|---|---|
| [demo_enet_cv_01.R](demo_enet_cv_01.R) | 10-fold CV on synthetic sparse-signal data for alpha in {0, 0.5, 1} (n=300, p=100), plus the ridge-only `lambda_2_lars` conversion scenario (n=300, p=200); each fit prints the CV-MSE curve (descending, glmnet-ordered lambda grid) and checks `lambda.min <= lambda.1se` and the curve/grid length agreement via `stopifnot()` |

### ElasticNetCV APIs used

- `ElasticNetCV$new()` with
  `$fit(X, y, alpha, n_folds, n_lambda, lambda_min_ratio, seed, standardize,
  intercept, max_iter, tol)`, `$cv_min()`, `$cv_1se()`, `$index_min()`,
  `$index_1se()`, `$get_lambdas()`, `$get_cv_errors()`, `$get_cv_std()`.

The C++ counterpart (`demo_mlm_ms_02_enet_cv_ccd`) additionally has a Part A
that fits standalone elastic-net regularization paths via `enet_gaussian`;
the R bindings expose only the CV class (no path class), so that part is
omitted here. The Python port
(`Python/ml_methods/model_selection/demo_enet_cv_01.py`) covers both parts.

### Running the ElasticNetCV demo

```bash
Rscript R/ml_methods/model_selection/demo_enet_cv_01.R
```

Console output only.

---

## Elastic-Net glmnet Cross-Check

Reference generator for cross-checking the C++ coordinate-descent elastic-net
solver — `trex::ml_methods::model_selection::elastic_net_gaussian` (path) and
`elastic_net_cv_gaussian` (K-fold CV) — against R's `glmnet` / `cv.glmnet`,
on the same sparse factor-model design used for the SPCA / tsolver
evaluations.

Why this is the strongest evidence for the lambda_2 machinery: `glmnet` *is*
cyclic coordinate descent, so the C++ CD engine should match glmnet's
coefficient path to solver tolerance (not merely within CV noise), unlike the
SVD-based ridge CV variants which only approximate `cv.glmnet`'s
`lambda.min` / `lambda.1se`.

### Files

| File | Description |
|---|---|
| [demo_en_glmnet_compare.R](demo_en_glmnet_compare.R) | Generates the factor-model design, runs `glmnet` / `cv.glmnet` for alpha in {0.0, 0.5, 1.0}, and dumps full-precision reference CSVs |
| `rdump_en/` | The reference CSV dumps consumed by the C++ comparison program |

### `rdump_en/` contents

| File(s) | Content |
|---|---|
| `meta.csv` | Single row: `n`, `p` |
| `Xn.csv`, `y.csv` | Centered factor design (n x p) and PC1 score (n x 1), full precision |
| `glmnet_lambda_<tag>.csv` | `fit$lambda` (descending), one column, per alpha tag `000` / `050` / `100` |
| `glmnet_beta_<tag>.csv` | p x nlambda coefficient path (original scale) |
| `glmnet_a0_<tag>.csv` | 1 x nlambda intercepts |
| `cv_glmnet_<tag>.csv` | `lambda.min`, `lambda.1se` (mean over CV repetitions), one row |

### Running the glmnet cross-check

Requires the CRAN `glmnet` package.

```bash
Rscript R/ml_methods/model_selection/demo_en_glmnet_compare.R
```

The dumps are written into `rdump_en/` next to the script; the C++ side then
reads them and compares its own path and CV results.

---

## Counterparts

- C++: [cpp/ml_methods/model_selection/](../../../cpp/ml_methods/model_selection/)
  — `demo_mlm_ms_01_ridge_cv_svd` (ridge K-fold CV via SVD, mirrored by
  `demo_ridge_cv_01.R`) and `demo_mlm_ms_02_enet_cv_ccd` (cross-validated
  elastic net via CCD, CV part mirrored by `demo_enet_cv_01.R`).
  The script header names the companion as
  `demo_mlm_ms_03_en_glmnet_rcompare.cpp`; no file of that name currently
  exists in the cpp examples tree.

---

**Last updated**: 2026-07-06
