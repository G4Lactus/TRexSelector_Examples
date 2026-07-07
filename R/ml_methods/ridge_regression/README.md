# Ridge Regression — R Demos

## Purpose

Demonstrates the **RidgeSolver** R6 class of the **TRexSelectorNeo** package:
solving the L2-penalized least-squares problem

```
beta_hat = argmin ||y - X beta||_2^2 + lambda ||beta||_2^2
```

in both the low-dimensional (primal, n >= p) and high-dimensional (dual,
n < p) regimes, plus a sweep over the penalty parameter.

---

## Demos

| File | Description |
|---|---|
| [demo_ridge_01_in_memory.R](demo_ridge_01_in_memory.R) | Primal-path solve recovering a dense linear model (n=500, p=50), dual-path solve with sparse ground truth (n=80, p=300), and a lambda sweep over 1e-4 ... 1e3 showing monotone coefficient shrinkage and growing training residual |

The demo encodes its expected properties as `stopifnot()` checks (estimate
close to the truth at low noise; `||beta_hat||_2` non-increasing and the
residual non-decreasing in lambda).

### APIs used

- `RidgeSolver$new()` with `$solve(X, y, lambda)`, returning the coefficient
  vector `beta_hat` (length p).

---

## Running

```bash
Rscript R/ml_methods/ridge_regression/demo_ridge_01_in_memory.R
```

Console output only.

---

## Counterparts

- C++: [cpp/ml_methods/ridge_regression/](../../../cpp/ml_methods/ridge_regression/)
  — `demo_mlm_ridge_01` (primal/dual ridge solve and lambda sweep).

---

**Last updated**: 2026-07-06
