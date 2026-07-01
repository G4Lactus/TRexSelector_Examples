# Ridge Regression: Demonstration Suite

## Overview

This folder contains C++ examples for the **ridge regression** implementation in the `ml_methods` module.

Ridge regression adds an $\ell_2$-penalty to least squares in order to stabilize estimation and shrink coefficients toward zero. The demo in this folder illustrates ridge fitting in both low-dimensional and high-dimensional settings, together with the effect of varying the penalty parameter.

---

## Key concept

Given a design matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$, response vector $\boldsymbol{y} \in \mathbb{R}^n$, and penalty parameter $\lambda \ge 0$, ridge regression estimates

$$
\widehat{\boldsymbol{\beta}}_{\mathrm{ridge}}
=
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\left\{
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda \lVert \boldsymbol{\beta} \rVert_2^2
\right\}.
$$

---

## Start here

| Folder | Purpose |
|------|---------|
| [demo_mlm_ridge_01/](demo_mlm_ridge_01/README.md) | Demonstrates ridge regression in primal and high-dimensional regimes, plus a lambda sweep |

---

## Running

```bash
./build/debug/bin/ml_methods/ridge_regression/demo_mlm_ridge_01/demo_mlm_ridge_01
```

---

**Last updated**: 2026-07-01
