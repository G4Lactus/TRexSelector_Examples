# Demo: `SVDSolver` Usage

## Purpose

This demo illustrates basic usage of the `SVDSolver` class for computing a truncated singular value decomposition.

It studies three scenarios:

1. a direct SVD path in the regime $n \ge p$,
2. a wide-matrix Gram-path case where $p > 2n$,
3. an orthogonality check for the returned singular vectors.

The demo is intended as a numerical sanity check and interface demonstration rather than a performance benchmark.

---

## Mathematical setup

For a matrix $\mathbf{X} \in \mathbb{R}^{n \times p}$, the singular value decomposition is

$$
\mathbf{X} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\top.
$$

If only the top $M$ singular values are retained, the solver returns a truncated decomposition

$$
\mathbf{X}_M = \mathbf{U}_M \boldsymbol{\Sigma}_M \mathbf{V}_M^\top,
$$

where

- $\mathbf{U}_M \in \mathbb{R}^{n \times M}$,
- $\boldsymbol{\Sigma}_M \in \mathbb{R}^{M \times M}$,
- $\mathbf{V}_M \in \mathbb{R}^{p \times M}$.

This decomposition gives the best rank-$M$ approximation to $\mathbf{X}$ in the Frobenius norm sense.

---

## Subproblem 1: Direct path

The first part of the demo generates a random Gaussian matrix with

$$
n = 300, \qquad p = 80, \qquad M = 10,
$$

and computes a truncated SVD in the regime $n \ge p$.

The code then prints the dimensions of $\mathbf{U}$, $\mathbf{S}$, and $\mathbf{V}$, reports the relative reconstruction
error

$$
\frac{\|\mathbf{X}_M - \mathbf{X}\|_F}{\|\mathbf{X}\|_F},
$$

and lists the top $M$ singular values.

---

## Subproblem 2: Gram path

The second part studies a wide matrix with

$$
n = 100, \qquad p = 500, \qquad M = 8,
$$

so that $p > 2n$ and the source explicitly labels this as the Gram-path case.

This scenario is useful because wide matrices often benefit from alternative internal linear-algebra strategies. As in
the first scenario, the demo reports output dimensions and the relative Frobenius reconstruction error of the truncated
approximation.

---

## Subproblem 3: Orthogonality check

The third part tests whether the returned singular vectors are numerically orthonormal. It uses a random matrix with

$$
n = 200, \qquad p = 150, \qquad M = 20,
$$

computes the truncated decomposition, and evaluates

$$
\|\mathbf{U}^\top \mathbf{U} - \mathbf{I}\|_F
\qquad \text{and} \qquad
\|\mathbf{V}^\top \mathbf{V} - \mathbf{I}\|_F
$$

If both quantities are sufficiently small, the demo reports that $\mathbf{U}$ and $\mathbf{V}$ are orthonormal up to
numerical precision.

---

## Running

```bash
./build/debug/bin/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01
```

The console output is organized into three sections matching the three subproblems above.

---

## What to look for

When reading the console output, focus on these points:

- whether the dimensions of $\mathbf{U}$, $\mathbf{S}$, and $\mathbf{V}$ match the requested truncated rank $M$,
- whether the singular values are reported in descending order,
- whether the truncated reconstruction error is numerically reasonable,
- whether the orthogonality errors $\|\mathbf{U}^\top \mathbf{U} - \mathbf{I}\|_F$ and
  $\|\mathbf{V}^\top \mathbf{V} - \mathbf{I}\|_F$ are close to zero.

The same SVD workflow is also implemented in
[Python](../../../../Python/ml_methods/svd/demo_svd_01_in_memory.py) and
[R](../../../../R/ml_methods/svd/demo_svd_01_in_memory.R).

---

## Technical notes

- The source labels the first scenario as the “Direct Path” for $n \ge p$ and the second as the “Gram Path” for the
  wide-matrix regime $p > 2n$.
- The demo does not include a dedicated rank-deficient example; the three scenarios use full-rank random Gaussian
  matrices.
- The reconstruction reported by the demo is a **truncated** reconstruction using only the top $M$ singular components,
  so the error is not expected to be near machine precision unless $M$ is large enough to capture essentially all
  variation.

---

**Last updated**: 2026-07-08
