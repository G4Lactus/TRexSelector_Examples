# PCA Demos (Python)

## Purpose

Demonstrate principal component analysis via the `trex_selector_neo.ml_methods`
`PCA` class — fitting on a data matrix, inspecting scores/loadings/explained
variance, restoring the in-place preprocessed matrix, and transforming new
data into the learned component space.

APIs used:

- `PCA(center=True, normalize=True)` — configure centering and normalization
- `PCA.fit(X, M)` — preprocess `X` in place and return a `PCAResult` with the
  score matrix `Z` (n x M), loading matrix `V` (p x M), and
  `explained_variance` (fraction of total variance per retained component)
- `PCA.transform(X_new)` — project new observations onto the learned loadings
- `PCA.restore()` — undo the in-place preprocessing of the fitted matrix

All matrices must be float64 and Fortran-ordered (`np.asfortranarray`);
`fit()` additionally requires a writeable array because it preprocesses `X`
in place. The in-place, zero-copy behavior is deliberate: the class works on
an `Eigen::Map` of the caller's buffer — the same interface that serves both
in-memory and memory-mapped data — because the library targets data volumes
that can exceed RAM and therefore never takes hidden copies of `X`. Pass an
explicit copy (`X.copy(order="F")`) where the original values are needed
afterwards, or call `restore()` to undo the preprocessing.

---

## Demos

| Demo | Description |
|---|---|
| [demo_pca_01_in_memory.py](demo_pca_01_in_memory.py) | Part A: basic fit (500 x 200, M = 10) with explained-variance summary (asserts fractions are decreasing and sum to <= 1); Part B: restore round-trip (300 x 100, M = 5) with max/mean error check; Part C: transform new data (train 400 x 80, test 50 x 80, M = 8) and verify `transform()` on the untouched training data reproduces the fit scores. |

Data is synthetic standard-normal in every part. RNG seeds match the C++
demo, but numpy's `default_rng` stream differs from `std::mt19937`, so
numbers agree qualitatively, not bitwise.

---

## Running the Demos

```bash
python demo_pca_01_in_memory.py
```

Output is printed to the console (assumes a Python environment with
`trex_selector_neo` and `numpy` installed).

---

## Counterparts

- C++: [cpp/ml_methods/pca/demo_mlm_pca_01](../../../cpp/ml_methods/pca/demo_mlm_pca_01/) (this script is a direct mirror)

---

**Last updated**: 2026-07-06
