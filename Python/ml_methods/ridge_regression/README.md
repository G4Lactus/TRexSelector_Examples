# Ridge Regression Demos (Python)

## Purpose

Demonstrate L2-regularized least squares via the
`trex_selector_neo.ml_methods` `RidgeSolver` class — solving
`min_beta ||y - X beta||_2^2 + lambda * ||beta||_2^2` in both the primal
(n >= p) and dual (n < p) regimes, and illustrating the shrinkage effect of
the penalty across a lambda sweep.

APIs used:

- `RidgeSolver.solve(X, y, lambda_val)` — static method; returns the ridge
  coefficient vector

`X` must be float64 and Fortran-ordered (`np.asfortranarray`). Cross-validated
lambda selection (`RidgeCV`) is not exposed in the Python bindings — see the
C++ `ml_methods/model_selection` demos for that workflow.

---

## Demos

| Demo | Description |
|---|---|
| [demo_ridge_01_in_memory.py](demo_ridge_01_in_memory.py) | Part A: primal path (500 x 50, lambda = 1.0) recovering a dense ground truth (asserts the coefficient error is small); Part B: dual path (80 x 300, lambda = 10.0) with a sparse ground truth in the n < p regime; Part C: lambda sweep (200 x 40, lambda in 1e-4 ... 1e3) printing coefficient error, residual, and `\|beta_hat\|_2`, and asserting that the coefficient norm shrinks and the training residual grows monotonically in lambda. |

Data is generated from a known linear model `y = X beta_true + noise` with
standard-normal design. RNG seeds match the C++ demo, but numpy's
`default_rng` stream differs from `std::mt19937`, so numbers agree
qualitatively, not bitwise.

---

## Running the Demos

```bash
python demo_ridge_01_in_memory.py
```

Output is printed to the console (assumes a Python environment with
`trex_selector_neo` and `numpy` installed).

---

## Counterparts

- C++: [cpp/ml_methods/ridge_regression/demo_mlm_ridge_01](../../../cpp/ml_methods/ridge_regression/demo_mlm_ridge_01/) (this script is a direct mirror)

---

**Last updated**: 2026-07-06
