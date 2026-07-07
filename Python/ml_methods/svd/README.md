# SVD Demos (Python)

## Purpose

Demonstrate truncated singular value decomposition via the
`trex_selector_neo.ml_methods` `SVDSolver` class — computing the rank-M
factorization `U @ diag(S) @ V.T`, measuring reconstruction error, and
checking orthonormality of the singular vectors.

APIs used:

- `SVDSolver()` — solver object
- `SVDSolver.compute(X, M)` — rank-M truncated SVD; returns an `SVDResult`
  with `U` (n x M left singular vectors), `S` (M singular values, decreasing),
  and `V` (p x M right singular vectors)

All matrices must be float64 and Fortran-ordered (`np.asfortranarray`).

---

## Demos

| Demo | Description |
|---|---|
| [demo_svd_01_in_memory.py](demo_svd_01_in_memory.py) | Part A: direct / thin-Jacobi path on a tall matrix (300 x 80, M = 10) with relative reconstruction error, top singular values, and a full-rank (M = p) exact-reconstruction assert; Part B: Gram / kernel path on a wide matrix (100 x 500, M = 8) where p > 2*n routes through the n x n Gram matrix; Part C: orthogonality check (200 x 150, M = 20) asserting `\|U^T U - I\|_F` and `\|V^T V - I\|_F` are tiny. |

Data is synthetic standard-normal in every part. RNG seeds match the C++
demo, but numpy's `default_rng` stream differs from `std::mt19937`, so
numbers agree qualitatively, not bitwise.

---

## Running the Demos

```bash
python demo_svd_01_in_memory.py
```

Output is printed to the console (assumes a Python environment with
`trex_selector_neo` and `numpy` installed).

---

## Counterparts

- C++: [cpp/ml_methods/svd/demo_mlm_svd_01](../../../cpp/ml_methods/svd/demo_mlm_svd_01/) (this script is a direct mirror)

---

**Last updated**: 2026-07-06
