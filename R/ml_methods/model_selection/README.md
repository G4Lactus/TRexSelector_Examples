# Model Selection — R Demos

This folder contains two model-selection demos: a **RidgeCV demo** (K-fold
cross-validation for the ridge penalty via the TRexSelectorNeo `RidgeCV` R6
class) and an **ElasticNet / ElasticNetCV demo** (path fitting plus K-fold CV
for the elastic net via coordinate descent).

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

## Elastic-Net Path Fitting and K-Fold CV (ElasticNet / ElasticNetCV)

Demonstrates glmnet-style pathwise cyclic coordinate descent: Part A fits
standalone regularization paths (`ElasticNet`), Part B does K-fold CV
(`ElasticNetCV`) with lambda selection across the alpha range
(ridge / elastic net / lasso) and the pure-ridge
`lambda_2_lars = lambda.1se * p / 2` conversion used internally by TRexGVS.

| File | Description |
|---|---|
| [demo_enet_cv_01.R](demo_enet_cv_01.R) | Part A: auto-grid paths for alpha in {0, 0.5, 1}, low-dim (n=300, p=100) and high-dim (n=200, p=500), with dev-ratio monotonicity and `fit_grid()` reproduction checks. Part B: 10-fold CV for the same alphas (n=300, p=100) plus the ridge-only `lambda_2_lars` conversion (n=300, p=200), printing the CV-MSE curve and checking `lambda.min <= lambda.1se` via `stopifnot()` |

### ElasticNet / ElasticNetCV APIs used

- `ElasticNet$new()` with `$fit(X, y, alpha)` / `$fit_grid(X, y, lambdas, alpha)`,
  `$get_coef()`, `$get_lambdas()`, `$get_dev_ratio()`, `$converged()`.
- `ElasticNetCV$new()` with
  `$fit(X, y, alpha, n_folds, n_lambda, lambda_min_ratio, seed, standardize,
  intercept, max_iter, tol)`, `$cv_min()`, `$cv_1se()`, `$index_min()`,
  `$index_1se()`, `$get_lambdas()`, `$get_cv_errors()`, `$get_cv_std()`.

Both parts of the C++ counterpart (`demo_mlm_ms_02_enet_cv_ccd`) are covered:
Part A fits standalone regularization paths via the `ElasticNet` class
(`$fit()` / `$fit_grid()`), Part B the K-fold CV via `ElasticNetCV`. The Python
port (`Python/ml_methods/model_selection/demo_enet_cv_01.py`) matches.

### Running the ElasticNetCV demo

```bash
Rscript R/ml_methods/model_selection/demo_enet_cv_01.R
```

Console output only.

---

> The glmnet cross-check that used to sit here (`demo_en_glmnet_compare.R` +
> `rdump_en/`) is local-testing material, not a demo. It moved to the
> TRexSelector library test suite
> (`TRexSelector/cpp/tests/validation/ml_methods/model_selection/`).

---

## Counterparts

- C++: [cpp/ml_methods/model_selection/](../../../cpp/ml_methods/model_selection/)
  — `demo_mlm_ms_01_ridge_cv_svd` (ridge K-fold CV via SVD, mirrored by
  `demo_ridge_cv_01.R`) and `demo_mlm_ms_02_enet_cv_ccd` (elastic-net path fit
  and K-fold CV via CCD, both parts mirrored by `demo_enet_cv_01.R`).
  The script header names the companion as
  `demo_mlm_ms_03_en_glmnet_rcompare.cpp`; no file of that name currently
  exists in the cpp examples tree.

---

**Last updated**: 2026-07-06
