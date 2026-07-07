# PCA — R Demos

## Purpose

Demonstrates the **Principal Component Analysis (PCA)** R6 class of the
**TRexSelectorNeo** package: fitting a truncated PC decomposition on Gaussian
data, inspecting explained variance, verifying the in-place restore
round-trip, and projecting new data into a fitted PC space.

---

## Demos

| File | Description |
|---|---|
| [demo_pca_01_in_memory.R](demo_pca_01_in_memory.R) | Basic PCA fit (score/loading dimensions, explained-variance profile), restore round-trip accuracy after in-place preprocessing, and transformation of new test data via `$transform()` |

The demo encodes its expected properties as `stopifnot()` checks
(non-increasing explained variance, exact restore round-trip, `transform()`
reproducing the fit scores).

### APIs used

- `PCA$new(X, M, center = TRUE, normalize = TRUE)` with `$fit()`,
  `$transform(X_new)`, `$restore()`, `$get_scores()`, `$get_loadings()`,
  `$get_explained_variance()`, `$get_means()`, `$get_norms()`.

### In-place (zero-copy) design

`PCA$new(X, ...)` maps the memory of the R matrix `X` directly — the same
`Eigen::Map` interface that serves both in-memory and memory-mapped data.
This is deliberate: the library targets data volumes that can exceed RAM, so
it never takes hidden copies of `X`. `$fit()` therefore standardizes `X`
**in place** and `$restore()` reverses the preprocessing exactly. Where a
pristine copy is needed, take an explicit deep copy with `X_orig <- X + 0`
(a plain alias shares the same buffer).

---

## Running

```bash
Rscript R/ml_methods/pca/demo_pca_01_in_memory.R
```

Console output only.

---

## Counterparts

- C++: [cpp/ml_methods/pca/](../../../cpp/ml_methods/pca/)
  — `demo_mlm_pca_01` (PCA fit, explained variance, restore round-trip,
  transform of new data).

---

**Last updated**: 2026-07-06
