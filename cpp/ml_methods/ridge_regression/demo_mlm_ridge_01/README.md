# Demo: Ridge Regression via `RidgeSolver`

## Purpose

This demo illustrates basic usage of the `RidgeSolver` class for $\ell_2$-penalized linear regression, and
studies three complementary scenarios:

1. a **primal-style** low-dimensional setting with $n \ge p$,
2. a **high-dimensional** setting with $n < p$,
3. a **lambda sweep** showing how the solution changes as the ridge penalty grows.

The demo is intended as an introductory numerical example rather than a full model-selection workflow.

---

## Statistical model and notation

Let

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the design matrix,
- $\boldsymbol{y} \in \mathbb{R}^n$ the response vector,
- $\boldsymbol{\beta} \in \mathbb{R}^p$ the regression coefficients,
- $\lambda \ge 0$ the ridge penalty.

The synthetic data are generated from a linear model of the form

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta}_{\mathrm{true}} + \boldsymbol{\varepsilon},
$$

with Gaussian predictors and additive Gaussian noise.

---

## Ridge objective

For a fixed penalty value $\lambda$, ridge regression solves

$$
\widehat{\boldsymbol{\beta}}(\lambda) =
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\left\lbrace
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda \lVert \boldsymbol{\beta} \rVert_2^2
\right\rbrace.
$$

When $\boldsymbol{X}^\top \boldsymbol{X}$ is used directly, this corresponds to the familiar closed-form expression

$$
\widehat{\boldsymbol{\beta}}_{\mathrm{ridge}} =
(\boldsymbol{X}^\top \boldsymbol{X} + \lambda \boldsymbol{I})^{-1}\boldsymbol{X}^\top \boldsymbol{y}.
$$

The effect of $\lambda$ is to shrink the coefficients toward zero while improving numerical stability relative to
unregularized least squares.

---

## Subproblem 1: Primal regime

The first scenario uses

$$
n = 500, \qquad p = 50, \qquad \lambda = 1.0,
$$

so the regression problem lies in the $n \ge p$ regime.

The code generates a dense Gaussian design, sets all true coefficients to 1, adds small Gaussian noise, solves the ridge
problem, and reports the coefficient error
$\lVert \widehat{\boldsymbol{\beta}} - \boldsymbol{\beta}_{\mathrm{true}} \rVert_2$ together with the residual norm
$\lVert \boldsymbol{y} - \boldsymbol{X}\widehat{\boldsymbol{\beta}} \rVert_2$.

---

## Subproblem 2: High-dimensional regime

The second scenario uses

$$
n = 80, \qquad p = 300, \qquad \lambda = 10.0,
$$

so the problem lies in the $n < p$ regime.

Here the true coefficient vector is sparse, with only the first 10 entries nonzero, and the demo again reports
coefficient error and residual norm after solving the ridge problem. This part is useful because it shows that the same
solver interface is used in a setting where regularization is especially important.

---

## Subproblem 3: Lambda sweep

The third scenario fixes a Gaussian regression problem with

$$
n = 200, \qquad p = 40,
$$

and then evaluates the solver over the penalty grid

$$
\lambda \in \lbrace10^{-4}, 10^{-2}, 0.1, 1, 10, 100, 1000\rbrace
$$

For each penalty value, the demo prints both the coefficient error and the residual norm, which lets the reader inspect
the trade-off between shrinkage and data fit as $\lambda$ increases.

---

## Running

```bash
./build/debug/bin/ml_methods/ridge_regression/demo_mlm_ridge_01/demo_mlm_ridge_01
```

The console output is organized into three sections matching the three scenarios above.

---

## What to look for

When reading the output, focus on these points:

- how the reported coefficient error differs between the $n \ge p$ and $n < p$ settings,
- how the residual norm changes across scenarios,
- how increasing $\lambda$ changes the balance between coefficient shrinkage and fit to the observed response,
- how the solver behaves under both dense and sparse ground-truth coefficient patterns.

The same ridge regression workflow is also implemented in
[Python](../../../../Python/ml_methods/ridge_regression/demo_ridge_01_in_memory.py) and
[R](../../../../R/ml_methods/ridge_regression/demo_ridge_01_in_memory.R).

---

## Technical notes

- The current source labels the first two cases as “Primal Path” and “Dual Path,” corresponding to the low-dimensional
  and high-dimensional regimes used in the demo output.
- The demo does not compute or print an OLS ($\lambda = 0$) baseline; it reports coefficient error and residual norm for
  the penalized solutions only.
- The examples use random Gaussian designs and report solver outputs rather than a dedicated ill-conditioning or
  collinearity study.

---

**Last updated**: 2026-07-09
