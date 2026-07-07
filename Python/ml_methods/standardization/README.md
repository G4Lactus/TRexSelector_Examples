# Standardization Demos (Python)

## Purpose

Demonstrate column standardization via the `trex_selector_neo.ml_methods`
scaler classes — fitting, in-place forward transformation, and in-place
inverse transformation (reconstruction check).

APIs used:

- `ZScoreScaler(with_mean=True, with_std=True)` — center columns and scale to
  unit variance
- `LpNormScaler(NormType.L2, True)` — center columns and scale to unit L2 norm
- `NormType` — norm selector for `LpNormScaler`
- Common methods: `fit()`, `transform_inplace()`, `inverse_transform_inplace()`
- `numpy_to_memmap` (from `trex_selector_neo.utils`, mmap variant only)

Both demos work on an 8 x 4 matrix drawn from N(3, 2^2) and print column-wise
mean, std (ddof=1), and L2 norm before/after scaling, plus the maximum
reconstruction error after the inverse transform.

---

## Demos

| Demo | Description |
|---|---|
| [demo_standardization_01_in_memory.py](demo_standardization_01_in_memory.py) | Part A: Z-scaling via `ZScoreScaler`; Part B: L2-norm scaling via `LpNormScaler`. `transform_inplace()` mutates the supplied array in place (and returns `None`) — the demo transforms a Fortran-ordered copy (`X.copy(order="F")`, mirroring R's `X + 0`) to preserve the original for the reconstruction check. |
| [demo_standardization_02_mmap.py](demo_standardization_02_mmap.py) | Same two scalers, but the matrix is first written to a temporary binary file via `numpy_to_memmap()`. The demo takes `X_mmap.to_numpy().copy(order="F")` (mirrors R's `as.matrix(X_mmap)`) because calling `transform_inplace()` on the zero-copy view directly would mutate the binary file on disk. |

---

## Running the Demos

```bash
python demo_standardization_01_in_memory.py
python demo_standardization_02_mmap.py
```

Output is printed to the console. The mmap demo cleans up its temporary
binary file on exit.

---

## Counterparts

- C++: [cpp/ml_methods/scaler_methods/demo_mlm_scaler_01](../../../cpp/ml_methods/scaler_methods/)
- R: `R/ml_methods/standardization/demo_standardization_01_in_memory.R` and
  `demo_standardization_02_mmap.R` (these scripts are direct mirrors)

---

**Last updated**: 2026-07-06
