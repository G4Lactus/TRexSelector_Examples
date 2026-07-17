# SVD: Demonstration Suite

## Overview

This folder contains C++ examples for the **singular value decomposition (SVD)** implementation in the `ml_methods`
module.

SVD is used in this project for stable matrix decompositions, dimensionality reduction, and regression-related
computations. In particular, it is relevant for PCA and for numerically stable ridge-regression workflows.

Equivalent SVD demos are also available in the Python and R packages; see
[Python/ml_methods/svd/](../../../Python/ml_methods/svd/README.md) and
[R/ml_methods/svd/](../../../R/ml_methods/svd/README.md) for the corresponding walkthroughs.

---

## Key concept

Given a matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$, the singular value decomposition is

$$
\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\top,
$$

where

- $\mathbf{U}$ contains left singular vectors,
- $\boldsymbol{\Sigma}$ contains singular values,
- $\mathbf{V}$ contains right singular vectors.

A truncated rank- $M$ approximation uses only the leading $M$ singular components:

$$
\mathbf{X}_M = \mathbf{U}_M \boldsymbol{\Sigma}_M \mathbf{V}_M^\top.
$$

---

## Start here

| Folder | Purpose |
| ------ | ------- |
| [demo_mlm_svd_01/](demo_mlm_svd_01/README.md) | Demonstrates basic `SVDSolver` usage, truncated reconstruction, and orthogonality checks |

---

## Running

```bash
./build/debug/bin/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01
```

---

**Last updated**: 2026-07-09
