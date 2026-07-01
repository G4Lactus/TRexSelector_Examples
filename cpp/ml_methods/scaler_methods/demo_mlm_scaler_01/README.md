# Demo: `ZScoreScaler` and `LpNormScaler`

## Purpose

This demo illustrates two column-scaling utilities used for data normalization:

- `ZScoreScaler`, which standardizes columns using a mean-and-scale transformation,
- `LpNormScaler`, which centers columns and rescales them to unit $L_p$-norm.

In addition to basic in-memory usage, the demo also covers inverse-transform accuracy, serialization of fitted scaler parameters, and scaling on a memory-mapped matrix.

---

## Mathematical setup

Let $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote a data matrix with columns $\boldsymbol{x}_1, \dots, \boldsymbol{x}_p$.

### Z-score scaling

For each column $j$, the transformed column is

$$
\tilde{\boldsymbol{x}}_j
=
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}{s_j},
$$

where $\bar{x}_j$ is the sample mean and $s_j$ is the stored scale parameter.

### Lp-norm scaling

With centering enabled, the transformed column is

$$
\tilde{\boldsymbol{x}}_j
=
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_p}.
$$
$$
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_p}.
$$

The demo explicitly uses both $p =2$ and $p =1$ in memory, and uses the $L_2$ version again in the memory-mapped example.

---

## Subproblem 1: In-memory scaler comparison

The first part of the demo generates a dense random matrix with

$$
n = 2000, \qquad p = 10000,
$$

fits three scalers on the same data, and then applies them to separate working copies:

- `ZScoreScaler`,
- `LpNormScaler` with $L_2$-norm,
- `LpNormScaler` with $L_1$-norm.

This section is mainly about showing how the scaler interfaces are fit once and then applied in place to large matrices.

---

## Subproblem 2: Inverse-transform accuracy

The second part tests whether scaling can be numerically reversed after an in-place transform. It generates a random matrix with

$$
n = 1000, \qquad p = 5000,
$$

then applies `fit()`, `transform_inplace()`, and `inverse_transform_inplace()` for both `ZScoreScaler` and the $L_2$-based `LpNormScaler`.

The source reports the maximum and mean absolute reconstruction errors after the round-trip, so this subproblem checks interface correctness and numerical reversibility rather than downstream modeling behavior.

---

## Subproblem 3: Serialization

The third part demonstrates persistence of a fitted scaler. A `ZScoreScaler` is fit on synthetic data with

$$
n = 1000, \qquad p = 500,
$$

saved to disk, loaded back into a fresh scaler object, and then compared through the norms of the differences in stored means and scales.

This section verifies that fitted preprocessing parameters can be serialized and restored without loss, which is useful when training and deployment happen in separate stages.

---

## Subproblem 4: Memory-mapped matrix

The fourth part applies scaling to a memory-mapped matrix with

$$
n = 20000, \qquad p = 100000.
$$

The source creates a writable mapped matrix, fills it deterministically, fits an $L_2$-based `LpNormScaler`, performs in-place transform and inverse transform operations, and then removes the backing file.

This subproblem shows that the scaler interface works with mapped matrix objects in addition to ordinary in-memory Eigen matrices.

---

## Running

```bash
./build/debug/bin/ml_methods/scaler_methods/demo_mlm_scaler_01/demo_mlm_scaler_01
```

The demo prints status messages for fitting, transforming, inverse transforming, serialization checks, and mapped-matrix cleanup.

---

## What to look for

When reading the console output, focus on these points:

- whether all scalers fit and transform successfully on the large in-memory matrix,
- whether inverse transforms return the data to the original values up to very small numerical error,
- whether serialized and reloaded `ZScoreScaler` parameters match essentially exactly,
- whether the memory-mapped example completes fit, transform, inverse transform, and cleanup without error.

---

## Technical notes

- The source uses `ZScoreScaler` together with `LpNormScaler` in both $L_1$ and $L_2$ modes in the in-memory comparison.
- The serialization example uses `ZScoreScaler::save()` and `ZScoreScaler::load()` and then compares stored means and scales numerically.
- The mapped-matrix example uses `LpNormScaler` with $L_2$-norm rather than `ZScoreScaler`.
- The fallback filenames in the current source are `scaler_test.bin` for serialization and `mmap_matrix.bin` for the mapped matrix, unless overridden by compile-time macros such as `TREX_DEMO_SCALER_SER_PATH` and `TREX_DEMO_MMAP_PATH`.

---

**Last updated**: 2026-07-01
