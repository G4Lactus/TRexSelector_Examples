# SVD — R Demos

## Purpose

Demonstrates the **SVDSolver** R6 class of the **TRexSelectorNeo** package:
truncated rank-M singular value decomposition `X_M = U_M diag(S_M) V_M^T` on
Gaussian data, covering both internal computational paths and the standard
numerical sanity checks.

---

## Demos

| File | Description |
|---|---|
| [demo_svd_01_in_memory.R](demo_svd_01_in_memory.R) | Truncated SVD on a tall matrix (direct path, n >= p) with full-rank reconstruction check, on a wide matrix (Gram/kernel path, p > 2n), and orthonormality checks `\|U^T U - I\|_F` / `\|V^T V - I\|_F` |

The demo encodes its expected properties as `stopifnot()` checks
(non-increasing singular values, near-zero full-rank reconstruction error,
orthonormal singular-vector blocks).

### APIs used

- `SVDSolver$new()` with `$compute(X, M)`, returning a list with elements
  `U` (n x M), `S` (length M), and `V` (p x M).

---

## Running

```bash
Rscript R/ml_methods/svd/demo_svd_01_in_memory.R
```

Console output only.

---

## Counterparts

- C++: [cpp/ml_methods/svd/](../../../cpp/ml_methods/svd/)
  — `demo_mlm_svd_01` (direct and Gram paths, truncated reconstruction,
  orthogonality checks).

---

**Last updated**: 2026-07-06
