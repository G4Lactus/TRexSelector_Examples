# Demo: In-Memory HAC Clustering Suite

## Purpose

This demo illustrates several in-memory uses of hierarchical agglomerative clustering (HAC) that are relevant for the TRexSelector project and, in particular, for grouped variable selection workflows.

We focus on six core subproblems, all of which are enabled and run in `main()`:

1. clustering samples with Euclidean distance,
2. clustering variables with correlation-based distance,
3. accelerating correlation-based clustering with Locality-Sensitive Hashing (LSH),
4. comparing linkage strategies on structured synthetic data (exact correlation),
5. comparing linkage strategies under the pure $O(1)$ LSH correlation approximation,
6. clustering AR(1) Toeplitz block structures, where correlation decays with distance.

Throughout the demo, clustering quality is assessed by the **Adjusted Rand Index (ARI)**, which compares the recovered cluster labels against the known ground-truth labels from the synthetic data generator.

---

## Statistical model and notation

Let

- $n$ the number of samples,
- $p$ the number of variables,
- $K$ the number of true clusters or blocks,
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $z_i \in \{1,\dots,K\}$ the true cluster label of sample $i$,
- $b(j) \in \{1,\dots,K\}$ the true block label of variable $j$.

The recovered clustering is evaluated with the Adjusted Rand Index

$$
\mathrm{ARI}
:=
\frac{\mathrm{Index} - \mathbb{E}[\mathrm{Index}]}
{\mathrm{MaxIndex} - \mathbb{E}[\mathrm{Index}]},
$$

where $\mathrm{ARI} = 1$ indicates perfect agreement between true and predicted labels.

---

## Subproblem 1: Sample clustering

### SP1 Data model

In the sample-clustering scenario, each row belongs to one of $K$ clusters and is generated according to

$$
X_{ij} = \mu_{z_i} + \sigma_{z_i} \varepsilon_{ij},
\qquad
\varepsilon_{ij} \sim \mathcal{N}(0,1).
$$

Here, $\mu_{z_i}$ is the cluster-specific mean and $\sigma_{z_i}$ is the cluster-specific scale.

### SP1 Mathematical approach

The goal is to cluster the rows $\boldsymbol{x}_1, \dots, \boldsymbol{x}_n$ using Euclidean distance,

$$
d(i,i')
=
\lVert \boldsymbol{x}_i - \boldsymbol{x}_{i'} \rVert_2.
$$

A **single-linkage** dendrogram is built with the SLINK algorithm. The dendrogram is then cut to recover $K$ clusters, and the predicted labels are compared with the true sample labels using ARI.

---

## Subproblem 2: Variable clustering with correlation distance

### SP2 Data model

In the variable-clustering scenario, variables are generated in correlated blocks.
Each variable $j$ belongs to a block $b(j)$, and the entries follow a latent-factor model of the form

$$
X_{ij}
=
\mu_{b(j)}
+
\sigma_{b(j)}
\left(
\alpha L_{i,b(j)} + \gamma \varepsilon_{ij}
\right),
$$

where

- $L_{i,b(j)} \sim \mathcal{N}(0,1)$ is the latent block signal,
- $\varepsilon_{ij} \sim \mathcal{N}(0,1)$ is independent noise,
- $\alpha$ controls the shared within-block signal,
- $\gamma$ controls the noise contribution.

This creates blockwise collinearity among variables.

### SP2 Standardization

Before correlation-based clustering, each column is standardized by

$$
\widetilde{\boldsymbol{x}}_j
=
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_2},
$$

so that variables are mean-centered and scaled to unit $\ell_2$-norm.

### SP2 Mathematical approach

After standardization, variables are clustered using the correlation distance

$$
d(j,\ell) = 1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
$$

Thus, variables with either strong positive or strong negative correlation are treated as close. The demo applies SLINK to this distance structure, builds a single-linkage dendrogram on the variables, and cuts the tree into $K$ groups.
The resulting block assignments are compared against the known true block labels.

---

## Subproblem 3: Variable clustering with LSH approximation

### SP3 Data model

The data model is the same block-correlated latent-factor model as in Subproblem 2, but the focus now shifts from exact clustering to **approximate similarity search** for large variable sets.

### SP3 Mathematical approach

This subproblem compares three correlation-based strategies:

1. **Exact correlation**, using the full correlation geometry.
2. **LSH filter**, where hashing acts as a gatekeeper before more exact comparisons.
3. **LSH approximation**, where correlation is approximated through a SimHash-style angular sketch.

Conceptually, the LSH approximation replaces expensive high-dimensional similarity comparisons by fast hash-based surrogates:

$$
\mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell)
\approx
\widehat{\mathrm{corr}}_{\mathrm{LSH}}(\boldsymbol{x}_j, \boldsymbol{x}_\ell).
$$

The purpose is to study whether the block structure is still recovered accurately when correlation distances are approximated rather than computed exactly.

---

## Subproblem 4: Linkage-method comparison

### SP4 Data model

This comparison is again carried out on structured synthetic data with known cluster labels, so that different linkage rules can be judged by their ARI and runtime.

### SP4 Mathematical approach

The main question is how clustering changes when the dendrogram update rule changes. The demo compares several linkage principles, including:

- **Single linkage**, based on nearest-neighbor inter-cluster distance,
- **Complete linkage**, based on farthest-neighbor inter-cluster distance,
- **Average linkage (UPGMA)**, based on average inter-cluster distance,
- **Weighted average linkage (WPGMA)**, based on weighted average inter-cluster distance,
- **Ward linkage**, based on minimum increase in within-cluster variance,
- and additional centroid/median-style variants, which don't have the reducibility property.

In abstract form, HAC iteratively merges two clusters $ A $ and $ B $ minimizing a linkage criterion

$$
d(A,B),
$$

where the definition of $d(A,B)$ depends on the chosen linkage rule. This subproblem studies how that choice affects both accuracy and computational behavior on the same synthetic structure. The comparison here is run under the **exact** Pearson-correlation distance, and it also exercises centroid- and median-style variants that lack the reducibility property.

---

## Subproblem 5: LSH-approximate linkage comparison

### SP5 Data model

This subproblem uses the same block-correlated latent-factor model as Subproblems 2–4, but on a larger, imbalanced configuration. Because the latent factors create dense, well-separated blocks in correlation space, the geometry is favorable for approximate hashing.

### SP5 Mathematical approach

The comparison repeats the linkage study of Subproblem 4, but replaces the exact correlation distance with the pure $O(1)$ SimHash-style LSH approximation

$$
d_{\mathrm{LSH}}(j,\ell)
\approx
1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
$$

Single (SLINK), complete, average (UPGMA), and weighted-average (WPGMA) linkage are each run under this approximation. The purpose is to show that, when the block geometry is well separated, the approximate distance still recovers the correct structure at a substantially reduced per-distance cost.

---

## Subproblem 6: AR(1) Toeplitz block clustering

### SP6 Data model

Here the variables are generated with a **block-diagonal AR(1) Toeplitz** correlation structure. Within each block, the correlation between variables decays with their index separation,

$$
\mathrm{corr}(x_{j}, x_{j+m}) \approx \rho^{\,m},
$$

so that a variable is strongly correlated with its immediate neighbors but nearly uncorrelated with distant variables in the same block.

### SP6 Mathematical approach

This chaining geometry is used to contrast linkage rules under the **exact** Pearson-correlation distance. Single linkage (SLINK) follows the correlation gradient and chains each block together, whereas complete and average linkage tend to fragment the gradually decaying blocks. LSH is deliberately not used here, because the coarse angular quantization destroys the continuous correlation gradient.

---

## Running the demo

```bash
./build/debug/bin/ml_methods/hac_clustering/demo_mlm_hac_01/demo_mlm_hac_01
```

The demo prints its main diagnostics to the console.

---

## What to look for

When reading the output, focus on the following points:

- whether the recovered clusters match the known synthetic structure,
- how strongly ARI depends on the data geometry,
- whether the LSH approximation preserves the dominant block structure,
- how the linkage rule changes both runtime and clustering agreement.

You can also compare selected results with the validation material in [validation_mlm_hac_01_rcompare](../../validation/hac_clustering/validation_mlm_hac_01_rcompare/README.md).

---

## Technical notes

- All six subproblems are enabled in `main()`; each is guarded by its own `if(true)` toggle so an individual part can be switched off while iterating on another.
- Subproblems 4 and 6 use the exact Pearson-correlation distance, while Subproblems 3 and 5 exercise the LSH filter and/or pure $O(1)$ LSH approximation.
- The synthetic dimensions vary per subproblem (for example, $P$ ranges up to $100{,}000$ variables), so runtimes and memory footprints differ across the six parts.

---

**Last updated**: 2026-07-08
