# Model Selection Demos (Python)

## Purpose

Demonstrate K-fold cross-validation for penalized regression via the
`trex_selector_neo.ml_methods` model-selection classes — lambda-grid
construction, the CV-MSE curve, the `lambda.min` / `lambda.1se` selections,
and (for the elastic net) coordinate-descent path fitting.

APIs used:

- `RidgeCV()` — SVD-based ridge K-fold CV with
  `fit(X, y, n_folds, n_lambda, lambda_ratio, seed)`
- `ElasticNet()` — glmnet-style CCD path solver with `fit()` (auto grid),
  `fit_grid()` (explicit grid), `coef()`, `intercepts()`, `lambdas()`,
  `dev_ratio()`, `converged()`, `n_nonconverged()`
- `ElasticNetCV()` — K-fold CV over a glmnet-style elastic-net grid with
  `fit(X, y, alpha, n_folds, n_lambda, lambda_min_ratio, seed, ...)`
- Common CV getters: `cv_min()`, `cv_1se()`, `index_min()`, `index_1se()`,
  `lambdas()`, `cv_mse()`, `cv_sem()`

All classes require float64, Fortran-ordered (column-major) arrays;
`np.asfortranarray()` provides the required layout. Both demos work on
synthetic sparse-signal data (`y = X beta_true + noise` with 10 active
coefficients) and use the same n / p / fold / seed values as their C++
counterparts; NumPy's PCG64 stream differs from `std::mt19937`, so results
match qualitatively, not bitwise.

---

## Demos

| Demo | Description |
|---|---|
| [demo_ridge_cv_01.py](demo_ridge_cv_01.py) | K-fold ridge CV via `RidgeCV` across three dimensional regimes — low-dimensional (n=1000, p=500, 10-fold), high-dimensional (n=300, p=500, 10-fold), and very high-dimensional (n=300, p=1000, 5-fold). Each scenario prints the CV-MSE curve (ascending lambda grid) and asserts `lambda.min <= lambda.1se`. The C++ counterpart's 4th, memory-mapped scenario is omitted (the Python binding only accepts in-memory arrays). |
| [demo_enet_cv_01.py](demo_enet_cv_01.py) | Part A: `ElasticNet` auto-grid paths for alpha in {0, 0.5, 1} (n=300, p=100) and {0, 1} (n=200, p=500), printing nnz / L1-norm / dev.ratio along the path, asserting a non-decreasing deviance ratio, and checking that `fit_grid()` at the auto-generated grid reproduces the auto path exactly. Part B: `ElasticNetCV` 10-fold CV per alpha (n=300, p=100) plus the pure-ridge `lambda_2_lars = lambda.1se * p / 2` conversion used by TRexGVS (n=300, p=200). The CV grid follows glmnet's descending order and may hold fewer than `n_lambda` points due to glmnet's early termination. |

---

## Running the Demos

```bash
python demo_ridge_cv_01.py
python demo_enet_cv_01.py
```

Output is printed to the console; the built-in `assert` checks make the
scripts fail loudly if the CV selections are inconsistent.

---

## Counterparts

- C++: [cpp/ml_methods/model_selection/](../../../cpp/ml_methods/model_selection/)
  — `demo_mlm_ms_01_ridge_cv_svd` (mirrored by `demo_ridge_cv_01.py`) and
  `demo_mlm_ms_02_enet_cv_ccd` (mirrored by `demo_enet_cv_01.py`)
- R: `R/ml_methods/model_selection/demo_ridge_cv_01.R` and
  `demo_enet_cv_01.R` (the R bindings expose no standalone elastic-net path
  class, so the R enet demo covers the CV part only)

---

**Last updated**: 2026-07-06
