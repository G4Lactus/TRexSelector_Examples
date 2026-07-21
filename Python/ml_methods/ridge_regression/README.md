# Ridge Regression: Demonstration Suite

## Overview

This folder contains Python examples for the **ridge regression** implementation in the
`trex_selector_neo.ml_methods` module.

Ridge regression adds an $\ell_2$-penalty to least squares in order to stabilize estimation and shrink coefficients
toward zero. The demo in this folder illustrates ridge fitting in both low-dimensional and high-dimensional settings,
together with the effect of varying the penalty parameter.

The C++ counterparts live in [cpp/ml_methods/ridge_regression/](../../../cpp/ml_methods/ridge_regression/); an
equivalent R walkthrough is available under `R/ml_methods/ridge_regression/`.

---

## Key concept

Given a design matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$, response vector $\boldsymbol{y} \in \mathbb{R}^n$,
and penalty parameter $\lambda \ge 0$, ridge regression estimates

$$
\widehat{\boldsymbol{\beta}}_{\mathrm{ridge}} =
\arg\min_{\boldsymbol{\beta} \in \mathbb{R}^p}
\left
\lVert \boldsymbol{y} - \boldsymbol{X}\boldsymbol{\beta} \rVert_2^2
+
\lambda \lVert \boldsymbol{\beta} \rVert_2^2
\right.
$$

---

## Start here

| Folder | Purpose |
| ------ | --------- |
| [demo_mlm_ridge_01/](demo_mlm_ridge_01/README.md) | Demonstrates ridge regression in primal and high-dimensional regimes, plus a lambda sweep |

---

## Running

```bash
python demo_mlm_ridge_01/demo_mlm_ridge_01.py
```

---

**Last updated**: 2026-07-21
