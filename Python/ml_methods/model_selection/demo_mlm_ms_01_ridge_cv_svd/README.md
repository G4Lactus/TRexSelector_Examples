# Demo: Ridge Regression K-Fold CV via the SVD-Based `RidgeCV`

## Purpose

This demo illustrates how to use `RidgeCV` for K-fold cross-validation in ridge regression across several synthetic
Gaussian regression settings, ranging from low-dimensional to high-dimensional regimes, and also includes a
memory-mapped design-matrix example. The focus is on selecting the regularization parameter $\lambda$ from
cross-validated prediction error while using a numerically stable SVD-based solver.

It mirrors the corresponding C++ example in
[cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/](../../../../cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/).

The demo is organized into four scenarios:

1. **Low-dimensional**: $n = 1000$, $p = 500$, 10-fold CV.
2. **High-dimensional**: $n = 300$, $p = 500$, 10-fold CV.
3. **High-dimensional**: $n = 300$, $p = 1000$, 5-fold CV.
4. **Memory-mapped dataset**: $n = 2000$, $p = 500$, 10-fold CV with the design matrix stored on disk
   (the C++ demo uses $n = 5000$; scaled down here to keep the Python run fast).

---

## Statistical model and notation

Let

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the design matrix,
- $\boldsymbol{y} \in \mathbb{R}^n$ the response vector,
- $\boldsymbol{\beta} \in \mathbb{R}^p$ the regression coefficients,
- $\lambda \ge 0$ the ridge penalty parameter.

The synthetic data follow the Gaussian linear model

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta}_{\mathrm{true}} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n).
$$

In the demo, the true coefficient vector is sparse: a fixed number of entries of $\boldsymbol{\beta}_{\mathrm{true}}$
are set to 1, placed at evenly spaced coordinates, while the remaining coefficients are zero.

---

## Ridge regression objective

For a given penalty parameter $\lambda$, ridge regression estimates the coefficients by solving

```math
\widehat{\boldsymbol{\beta}}(\lambda)
=
\underset{\boldsymbol{\beta} \in \mathbb{R}^p}{\arg\min}
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+ \lambda \lVert \boldsymbol{\beta} \rVert_2^2
.
```

In the implementation, the CV routine performs per-fold centering and column $\ell_2$-normalization and uses an SVD
(`Eigen::JacobiSVD` on the C++ side) for the foldwise numerical solve.

---

## Cross-validation criterion

The penalty parameter is selected by K-fold cross-validation. For each candidate value of $\lambda$, the data are
partitioned into $K$ folds, the model is fit on $K-1$ folds, and prediction error is evaluated on the held-out fold.

This yields a cross-validated mean squared error curve

$$
\mathrm{CV\text{-}MSE}(\lambda),
$$

from which the demo reports two common choices:

- $\lambda_{\min}$: the value minimizing cross-validated MSE,
- $\lambda_{1\mathrm{se}}$: the largest value of $\lambda$ within one standard error (SE) of the minimum.

The script explicitly prints both `lambda.min` and `lambda.1se` and checks the expected ordering `lambda.min <=
lambda.1se` in each scenario.

---

## Scenario 1: Low-dimensional ridge CV

In the first scenario, the design satisfies $n > p$, so the regression problem is overdetermined. Mathematically, this
is the most classical case for ridge regression: the penalty is mainly used for shrinkage and stability rather than for
handling rank deficiency. The demo reports the CV-MSE curve and the selected values of $\lambda_{\min}$ and
$\lambda_{1\mathrm{se}}$.

---

## Scenario 2: High-dimensional ridge CV

In the second scenario, the design satisfies $p > n$, so the regression problem becomes high-dimensional. In this
regime, regularization is more central because the unpenalized least-squares problem is no longer well posed in the
usual full-rank sense.

The ridge penalty ensures a stable solution even when the predictor dimension exceeds the sample size. This scenario
shows how the CV-based choice of $\lambda$ adapts once the problem enters the $p > n$ regime.

---

## Scenario 3: Very high-dimensional ridge CV

In the third scenario, the design is even more high-dimensional, with $p \gg n$. The demo reduces the number of folds
from 10 to 5 so that the training partitions remain sufficiently large relative to the sample size.

From a mathematical perspective, this is still the same ridge objective, but the estimation problem is more strongly
underdetermined without regularization. The scenario is useful for checking whether the cross-validation workflow
remains stable when the feature dimension grows further.

---

## Scenario 4: Memory-mapped design matrix

In the fourth scenario, the design matrix is first written to a memory-mapped backing file and then re-opened for
fitting. The important idea is that `RidgeCV.fit()` accepts any Fortran-ordered `float64` array, so the zero-copy
`to_numpy()` view of a `MemoryMappedMatrix` binds without copying the full design matrix into a separate dense object
at the interface boundary.

Two adaptations relative to the C++ demo: the C++ version re-opens the file **read-only** through a `ConstMapType`,
whereas the Python zero-copy view requires `ReadWrite` mode; and the sample size is scaled from $n = 5000$ to
$n = 2000$ to keep the demo fast. Note that the example is not a fully out-of-core approach because individual
train/test submatrices are still materialized in RAM inside the CV routine.

---

## Running the demo

```bash
python demo_mlm_ms_01_ridge_cv_svd.py
```

The demo prints the cross-validated MSE curve, the selected `lambda.min` and `lambda.1se` values, and an ordering check
for each scenario.

---

## What to look for

When reading the console output, check the following points:

- how the CV-MSE curve changes across the four dimensionality regimes,
- whether the reported `lambda.min` and `lambda.1se` values are sensible and satisfy the expected ordering,
- how strongly the selected ridge penalty changes when moving from $n > p$ to $p > n$ and then to $p \gg n$,
- how the memory-mapped scenario behaves relative to the in-memory cases.

In the memory-mapped scenario, also check that the temporary backing file is created, re-opened for fitting, and
removed successfully at the end.

---

## Technical notes

- `RidgeCV.fit(X, y, n_folds, n_lambda, lambda_ratio, seed)` requires `float64`, Fortran-ordered arrays
  (`np.asfortranarray`); the demo uses a 100-point lambda grid spanning a 1000x range, mirroring the C++ fit.
- The lambda grid is anchored internally using a scale derived from the standardized design and response.
- The memory-mapped backing file is created with `tempfile.mkstemp()` in the system temporary directory (the C++ demo
  uses the fallback filename `trex_mlm_demo_ms01_svd_mmap.bin`) and is removed in a `try/finally` block.
- The scratch file is only relevant to the fourth scenario.
- Data seeds match the C++ demo, but NumPy's PCG64 stream differs from `std::mt19937`, so results match qualitatively,
  not bitwise.

---

**Last updated**: 2026-07-21
