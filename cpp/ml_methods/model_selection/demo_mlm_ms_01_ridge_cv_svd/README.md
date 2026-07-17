# Demo: Ridge Regression K-Fold CV via JacobiSVD

## Purpose

This demo illustrates how to use `ridge_cv_svd` for K-fold cross-validation in ridge regression across several synthetic
Gaussian regression settings, ranging from low-dimensional to high-dimensional regimes, and also includes a
memory-mapped design-matrix example. The focus is on selecting the regularization parameter $\lambda$ from
cross-validated prediction error while using a numerically stable SVD-based solver.

The demo is organized into four scenarios:

1. **Low-dimensional**: $n = 1000$, $p = 500$, 10-fold CV.
2. **High-dimensional**: $n = 300$, $p = 500$, 10-fold CV.
3. **High-dimensional**: $n = 300$, $p = 1000$, 5-fold CV.
4. **Memory-mapped dataset**: $n = 5000$, $p = 500$, 10-fold CV with the design matrix stored on disk.

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

$$
\widehat{\boldsymbol{\beta}}(\lambda) =
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\left\lbrace\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2 +
\lambda \lVert \boldsymbol{\beta} \rVert_2^2 \right\rbrace.
$$

In the implementation, `ridge_cv_svd` performs per-fold centering and column $\ell_2$-normalization and uses
`Eigen::JacobiSVD` for the foldwise numerical solve.

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

The source explicitly prints both `lambda.min` and `lambda.1se` and checks the expected ordering `lambda.min <=
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

In the fourth scenario, the design matrix is first written to a memory-mapped backing file and then reopened in
read-only mode for fitting. The important idea is that `ridge_cv_svd::fit()` accepts an Eigen-compatible reference type,
so the mapped matrix can be consumed without copying the full design matrix into a separate dense object at the
interface boundary.

Mathematically, the regression and cross-validation problem are unchanged:

$$
\widehat{\boldsymbol{\beta}}(\lambda) =
\arg\min_{\boldsymbol{\beta}}
\left\lbrace
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda \lVert \boldsymbol{\beta} \rVert_2^2
\right\rbrace.
$$

What changes is the storage layer. The source notes, however, that this is not a fully out-of-core foldwise solver,
because individual train/test submatrices are still materialized in RAM inside the CV routine.

---

## Running the demo

```bash
./build/debug/bin/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/demo_mlm_ms_01_ridge_cv_svd
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

In the memory-mapped scenario, also check that the temporary backing file is created, reopened read-only for fitting,
and removed successfully at the end.

The same ridge CV workflow is also implemented in
[Python](../../../../Python/ml_methods/model_selection/demo_ridge_cv_01.py) and
[R](../../../../R/ml_methods/model_selection/demo_ridge_cv_01.R).

---

## Technical notes

- `ridge_cv_svd` uses `Eigen::JacobiSVD` rather than `BDCSVD` for numerical safety in the implementation shown here.
- The lambda grid is anchored internally using a scale derived from the standardized design and response, as described
  in the source comments.
- The memory-mapped file path is injected through `TREX_DEMO_MMAP_PATH` when provided, with a relative fallback filename
  otherwise.
- The scratch file is only relevant to the fourth scenario.

---

**Last updated**: 2026-07-09
