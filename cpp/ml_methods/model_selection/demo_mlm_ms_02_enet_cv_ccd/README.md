# Demo: Elastic-Net Path Fitting and K-Fold CV via Coordinate Descent

## Purpose

This demo illustrates how to use `enet_gaussian` for elastic-net path fitting and `enet_cv_ccd` for K-fold
cross-validation.

It studies Gaussian regression with elastic-net regularization across several synthetic-data settings and shows how the
regularization path changes with the elastic-net mixing parameter $\alpha$. In addition, it demonstrates how
cross-validation selects $\lambda$ for ridge, elastic-net, and lasso variants, including the ridge-specific conversion
used internally by TRexGVS.

The demo is organized into four scenarios:

1. **Low-dimensional path fitting**: regularization paths for $\alpha \in \{0, 0.5, 1\}$ with $n = 300$, $p = 100$.
2. **High-dimensional path fitting**: regularization paths for $\alpha \in \{0, 1\}$ with $n = 200$, $p = 500$.
3. **Low-dimensional cross-validation**: 10-fold CV for $\alpha \in \{0, 0.5, 1\}$ with $n = 300$, $p = 100$.
4. **Ridge CV with TRexGVS conversion**: 10-fold CV for $\alpha = 0$ with $n = 300$, $p = 200$, and explicit conversion
   from glmnet-scale $\lambda$ to the `lambda_2_lars` scale used downstream.

---

## Statistical model and notation

Let

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the design matrix,
- $\boldsymbol{y} \in \mathbb{R}^n$ the response vector,
- $\boldsymbol{\beta} \in \mathbb{R}^p$ the regression coefficients,
- $\lambda \ge 0$ the penalty parameter,
- $\alpha \in [0,1]$ the elastic-net mixing parameter.

The synthetic data follow the Gaussian linear model

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta}_{\mathrm{true}} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n).
$$

As in the ridge demo, the true signal is sparse: a fixed number of coefficients in $\boldsymbol{\beta}_{\mathrm{true}}$
are set to 1 at evenly spaced positions, while the remaining entries are zero.

---

## Elastic-net objective

For fixed $\alpha$ and $\lambda$, the elastic-net estimator solves

$$
\widehat{\boldsymbol{\beta}}(\lambda,\alpha)
=
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\left\{
\frac{1}{2n}\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda
\left[
\alpha \lVert \boldsymbol{\beta} \rVert_1
+
\frac{1-\alpha}{2}\lVert \boldsymbol{\beta} \rVert_2^2
\right]
\right\}.
$$

Here the residual sum of squares carries glmnet's $1/(2n)$ normalization (the
implementation uses `invN = 1/n`), which is the scale on which the demo's lambda
grid and the `lambda_2_lars = lambda * p / 2` conversion are defined.

The parameter $\alpha$ controls the penalty type:

- $\alpha = 0$: ridge regression,
- $\alpha = 1$: lasso,
- $0 < \alpha < 1$: elastic net.

The source implements pathwise cyclic coordinate descent and states that `enet_gaussian` follows glmnet-style semantics,
including a sequential strong rule and early termination logic.

---

## Path fitting

For a fixed value of $\alpha$, `enet_gaussian` computes a decreasing path of penalty values

$$
\lambda_1 > \lambda_2 > \cdots > \lambda_K,
$$

and returns the corresponding sequence of fitted coefficient vectors

$$
\widehat{\boldsymbol{\beta}}(\lambda_1,\alpha),
\widehat{\boldsymbol{\beta}}(\lambda_2,\alpha),
\dots,
\widehat{\boldsymbol{\beta}}(\lambda_K,\alpha).
$$

The demo summarizes each path by reporting selected $ \lambda $ values together with:

- the number of nonzero coefficients,
- the $ \ell_1 $-norm of the fitted coefficient vector,
- the deviance ratio.

This makes it possible to inspect how sparsity and fit quality evolve along the regularization path.

---

## Cross-validation criterion

For a fixed $\alpha$, `enet_cv_ccd` uses K-fold cross-validation to select the penalty parameter $\lambda$. The source
notes that the full-data lambda grid is constructed once and then reused across folds, following glmnet-style semantics.

As in the ridge demo, the main outputs are:

- $\lambda_{\min}$: the value minimizing cross-validated MSE,
- $\lambda_{1\mathrm{se}}$: the largest value within one standard error of the minimum.

The demo prints both quantities, the associated CV-MSE values, and an ordering check confirming that
$\lambda_{\min} \le \lambda_{1\mathrm{se}}$.

---

## Scenario 1: Low-dimensional path fitting

The first scenario studies the regularization path in a low-dimensional regime with $n = 300$ and $p = 100$. The
demo compares three values of the mixing parameter:

$$
\alpha \in \{0, 0.5, 1\}.
$$

This setting is useful for seeing how ridge, intermediate elastic-net shrinkage, and lasso differ in coefficient
sparsity and path shape when the number of predictors is moderate relative to the sample size.

---

## Scenario 2: High-dimensional path fitting

The second scenario studies the path in a high-dimensional regime with $n = 200$ and $p = 500$, so $p > n$. Here
the demo focuses on ridge and lasso, for which

$$
\alpha \in \{0, 1\}.
$$

This subproblem highlights how the path behaves when the regression problem becomes high-dimensional and regularization
is essential for stable estimation. It is especially useful for contrasting dense ridge-type shrinkage with sparse
lasso-type shrinkage.

---

## Scenario 3: Low-dimensional cross-validation

The third scenario performs 10-fold cross-validation in the low-dimensional setting $n = 300$, $p = 100$ for

$$
\alpha \in \{0, 0.5, 1\}.
$$

The goal is to compare how the selected $\lambda$ values differ across ridge, elastic net, and lasso when the same
underlying signal and noise structure are used. The output includes the CV-MSE curve, `lambda.min`, `lambda.1se`, and a
consistency check on the expected ordering of the selected penalties.

---

## Scenario 4: Ridge CV and TRexGVS conversion

The fourth scenario focuses on the pure-ridge case $\alpha = 0$, which is especially relevant for downstream TRexGVS
workflows. The source explains that for ridge, the returned cross-validated penalty is on the glmnet-reported scale and
is converted internally via

$$
\lambda_{2,\mathrm{lars}} = \lambda_{1\mathrm{se}} \cdot \frac{p}{2}.
$$

The demo prints both the `lambda.min` and `lambda.1se` values and then explicitly computes the converted ridge-scale
quantities

$$
\lambda_{2,\mathrm{lars}}^{(\min)} = \lambda_{\min} \cdot \frac{p}{2},
\qquad
\lambda_{2,\mathrm{lars}}^{(1\mathrm{se})} = \lambda_{1\mathrm{se}} \cdot \frac{p}{2},
$$

showing how the CV-selected penalty is mapped into the scale used by TRexGVS.

---

## Running the demo

```bash
./build/debug/bin/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/demo_mlm_ms_02_enet_cv_ccd
```

The demo prints condensed summaries of regularization paths and cross-validation results to the console.

---

## What to look for

When reading the console output, focus on the following points:

- how the number of nonzero coefficients changes along the path as $\lambda$ decreases,
- how the path differs between ridge, elastic net, and lasso,
- how the selected `lambda.min` and `lambda.1se` values vary across $\alpha$,
- how the high-dimensional regime changes the path behavior relative to the low-dimensional case,
- how the ridge-specific `lambda_2_lars = lambda * p / 2` conversion is reported in the final scenario.

The same elastic-net path fitting and CV workflow is also implemented in
[Python](../../../../Python/ml_methods/model_selection/demo_enet_cv_01.py) and
[R](../../../../R/ml_methods/model_selection/demo_enet_cv_01.R).

---

## Technical notes

- `enet_gaussian` uses pathwise cyclic coordinate descent with glmnet-style path construction and stopping behavior.
- `enet_cv_ccd` builds the lambda grid once on the full dataset and reuses it across folds, rather than rebuilding a
  separate grid for every fold.
- The demo is intended to explain path behavior and hyperparameter selection, not to benchmark prediction accuracy
  across many replicated simulations.

---

**Last updated**: 2026-07-09
