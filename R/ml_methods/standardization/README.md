# Standardization — R Demos

## Purpose

Demonstrates the column-standardization scalers of the **TRexSelectorNeo**
package on a small 8 x 4 Gaussian matrix, on both in-memory and memory-mapped
data.

---

## Demos

| File | Description |
|---|---|
| [demo_standardization_01_in_memory.R](demo_standardization_01_in_memory.R) | Z-scaling via `ZScoreScaler` (center + scale to unit variance) and L2-norm scaling via `LpNormScaler` (center + scale to unit L2 norm) on in-memory data, with forward and inverse transformations |
| [demo_standardization_02_mmap.R](demo_standardization_02_mmap.R) | Same scalers applied to mmap-backed data: the matrix is first written to disk via `convert_to_memory_mapped()`, then transformed in place |

Both demos print per-column mean / sd / L2-norm summaries before and after
each transformation and verify the inverse-transform round trip.

### APIs used

- `ZScoreScaler$new(center = TRUE, scale = TRUE)` with `$fit()`,
  `$transform_inplace()`, `$inverse_transform_inplace()`.
- `LpNormScaler$new(norm_type = 2L, center = TRUE)` with the same method set.
- `convert_to_memory_mapped()` (demo 02).

---

## Running

```bash
Rscript R/ml_methods/standardization/demo_standardization_01_in_memory.R
Rscript R/ml_methods/standardization/demo_standardization_02_mmap.R
```

Console output only.

---

## Counterparts

- C++: [cpp/ml_methods/scaler_methods/](../../../cpp/ml_methods/scaler_methods/)
  — `demo_mlm_scaler_01` (Z-score and Lp-norm column normalization).

---

**Last updated**: 2026-07-06
