# Demo: PCA Fit, Restore, and Transform

## Purpose

This demo illustrates how to use the PCA implementation in three basic workflows:

1. fitting PCA on a data matrix,
2. restoring the original matrix after in-place preprocessing,
3. transforming new data into the learned principal-component space.

The demo is designed as a practical introduction to the PCA interface and to the main matrix objects returned by a fit: principal-component scores, loading vectors, and explained variances.

---

## Statistical model and notation

Let

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $M \le p$ the number of retained principal components,
- $\boldsymbol{Z} \in \mathbb{R}^{n \times M}$ the score matrix,
- $\boldsymbol{V} \in \mathbb{R}^{p \times M}$ the loading matrix.

PCA constructs an orthogonal low-dimensional representation of the data by projecting the original variables onto $M$ principal directions.
In matrix form, the fitted scores are obtained by projecting the centered data onto the retained loading vectors.

---

## PCA representation

Conceptually, PCA approximates the data matrix by

$$
\boldsymbol{X} \approx \boldsymbol{Z}\boldsymbol{V}^\top,
$$

where the columns of $\boldsymbol{V}$ are the retained principal directions and the columns of $\boldsymbol{Z}$ are the corresponding principal-component scores.

The explained variances quantify how much of the retained variation is captured by each component. The demo reports the explained variance percentages for the retained $M$ components after fitting.

---

## Subproblem 1: Basic fit

The first part of the demo generates a random matrix with

$$
n = 500, \qquad p = 200, \qquad M = 10,
$$

fits PCA, and prints the dimensions of the returned score matrix $\boldsymbol{Z}$, loading matrix $\boldsymbol{V}$, and explained-variance vector $\boldsymbol{\lambda}$.

This part is mainly about interface understanding: after fitting, the user can inspect the shape of the PCA outputs and the relative contribution of each retained component to the retained variance budget.

---

## Subproblem 2: Restore round-trip

The second part studies whether the PCA object can safely undo its in-place preprocessing operations. The code generates a random matrix with

$$
n = 300, \qquad p = 100, \qquad M = 5,
$$

fits PCA on a working copy, restores the matrix, and then compares the restored matrix against the original one using maximum and mean absolute reconstruction errors.

The point here is not low-rank approximation error, but interface correctness: `restore()` is checked as a round-trip operation on the preprocessed training matrix representation.

---

## Subproblem 3: Transform new data

The third part separates training and test data. PCA is fit on a training matrix with

$$
n_{\mathrm{train}} = 400, \qquad p = 80, \qquad M = 8,
$$

and then applied to a test matrix with

$$
n_{\mathrm{test}} = 50
$$

to produce principal-component scores for unseen observations.

The demo also verifies that applying `transform()` to an untouched copy of the original training data reproduces the same score matrix returned during fitting, up to numerical tolerance.

---

## Running the demo

```bash
./build/debug/bin/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01
```

The demo prints matrix dimensions, explained-variance summaries, round-trip reconstruction errors, and transform consistency checks to the console.

---

## What to look for

When reading the console output, focus on the following points:

- whether the score matrix \( \boldsymbol{Z} \) and loading matrix \( \boldsymbol{V} \) have the expected dimensions,
- whether the explained-variance entries are reported consistently for the retained components,
- whether `restore()` returns the working matrix to its original state up to very small numerical error,
- whether `transform()` produces the expected output shape for new data,
- whether transforming the original training data reproduces the fitted scores.

---

## Technical notes

- The PCA object is constructed from an `Eigen::Map`, so the training matrix is passed through an Eigen-compatible interface.
- The demo uses synthetic Gaussian data generated from standard normal draws in each subproblem.
- The current README sentence “explained variance ratio per component (should sum to 1 across all components)” is only accurate when interpreted over the retained components printed by the demo, not over the full original $p$-dimensional variance decomposition.

---

**Last updated**: 2026-07-01
