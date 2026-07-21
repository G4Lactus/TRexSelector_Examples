# Demo: In-Memory HAC Clustering Suite

## Purpose

This demo illustrates several in-memory uses of hierarchical agglomerative clustering (HAC) that are relevant for the
TRexSelector project and, in particular, for grouped variable selection workflows.

We focus on six core subproblems, all of which are enabled and run in the main block:

1. clustering samples with Euclidean distance,
2. clustering variables with correlation-based distance,
3. accelerating correlation-based clustering with Locality-Sensitive Hashing (LSH),
4. comparing linkage strategies on structured synthetic data (exact correlation),
5. comparing linkage strategies under the pure $O(1)$ LSH correlation approximation,
6. clustering AR(1) Toeplitz block structures, where correlation decays with distance.

Throughout the demo, clustering quality is assessed by the **Adjusted Rand Index (ARI)**, which compares the recovered
cluster labels against the known ground-truth labels from the synthetic data generator.

It mirrors the corresponding C++ example in
[cpp/ml_methods/hac_clustering/demo_mlm_hac_01/](../../../../cpp/ml_methods/hac_clustering/demo_mlm_hac_01/). The C++
demo runs at large scale (up to $p = 100{,}000$ variables with OpenMP-parallel data generation); the Python mirror
keeps the same six scenarios, generators, and evaluation pipeline at sizes scaled down by roughly a factor of 100 so
the demo stays fast (the per-scenario scale factors are noted in the script header).

---

## Statistical model and notation

Let

- $n$ be the number of samples,
- $p$ be the number of variables,
- $K$ be the number of true clusters or blocks,
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $z_i \in \lbrace1,\dots,K\rbrace$ the true cluster label of sample $i$,
- $b(j) \in \lbrace1,\dots,K\rbrace$ the true block label of variable $j$.

The recovered clustering is evaluated with the Adjusted Rand Index

$$
\mathrm{ARI} :=
\frac{\mathrm{Index} - \mathbb{E}\lbrack\mathrm{Index}\rbrack}
{\mathrm{MaxIndex} - \mathbb{E}\lbrack\mathrm{Index}\rbrack},
$$

where $\mathrm{ARI} = 1$ indicates perfect agreement between true and predicted labels.

---

## Subproblem 1: Sample clustering

### SP1 Data model

In the sample-clustering scenario ($n = 300$, $p = 1000$, $K = 3$; cpp: $n = 1500$, $p = 50000$), each row belongs to
one of $K$ clusters and is generated according to

$$
X_{ij} = \mu_{z_i} + \sigma_{z_i} \varepsilon_{ij},
\qquad
\varepsilon_{ij} \sim \mathcal{N}(0,1).
$$

Here, $\mu_{z_i}$ is the cluster-specific mean and $\sigma_{z_i}$ is the cluster-specific scale.

### SP1 Mathematical approach

The goal is to cluster the rows $\boldsymbol{x}_1, \dots, \boldsymbol{x}_n$ using Euclidean distance,

$$
d(i,i') = \lVert \boldsymbol{x}_i - \boldsymbol{x}_{i'} \rVert_2 .
$$

A **single-linkage** dendrogram is built with the SLINK algorithm. The dendrogram is then cut to recover $K$ clusters,
and the predicted labels are compared with the true sample labels using ARI.

---

## Subproblem 2: Variable clustering with correlation distance

### SP2 Data model

In the variable-clustering scenario ($n = 500$, $p = 500$, $K = 5$; cpp: $n = 10000$, $p = 50000$), variables are
generated in correlated blocks.
Each variable $j$ belongs to a block $b(j)$, and the entries follow a latent-factor model of the form

$$
X_{ij} =
\mu_{b(j)} +
\sigma_{b(j)}
\left(
\alpha L_{i,b(j)} + \gamma \varepsilon_{ij}
\right)
$$

where

- $L_{i,b(j)} \sim \mathcal{N}(0,1)$ is the latent block signal,
- $\varepsilon_{ij} \sim \mathcal{N}(0,1)$ is independent noise,
- $\alpha = 0.85$ controls the shared within-block signal,
- $\gamma = 0.15$ controls the noise contribution.

This creates blockwise collinearity among variables.

### SP2 Standardization

Before correlation-based clustering, each column is standardized by

$$
\widetilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_2}
$$

so that variables are mean-centered and scaled to unit $\ell_2$-norm.

### SP2 Mathematical approach

After standardization, variables are clustered using the correlation distance

$$
d(j,\ell) = 1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
$$

Thus, variables with either strong positive or strong negative correlation are treated as close. The demo applies SLINK
to this distance structure, builds a single-linkage dendrogram on the variables, and cuts the tree into $K$ groups. The
resulting block assignments are compared against the known true block labels.

---

## Subproblem 3: Variable clustering with LSH approximation

### SP3 Data model

The data model is the same block-correlated latent-factor model as in Subproblem 2 ($n = 800$, $p = 500$; cpp:
$n = 40000$, $p = 50000$), but the focus now shifts from exact clustering to **approximate similarity search** for
large variable sets.

### SP3 Mathematical approach

This subproblem compares three correlation-based strategies:

1. **Exact correlation**, using the full correlation geometry.
2. **LSH filter**, where hashing acts as a gatekeeper before more exact comparisons.
3. **LSH approximation**, where correlation is approximated through a SimHash-style angular sketch.

Conceptually, the LSH approximation replaces expensive high-dimensional similarity comparisons by fast hash-based
surrogates:

$$
\mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell)
\approx
\widehat{\mathrm{corr}}_{\mathrm{LSH}}(\boldsymbol{x}_j, \boldsymbol{x}_\ell).
$$

The purpose is to study whether the block structure is still recovered accurately when correlation distances are
approximated rather than computed exactly.

---

## Subproblem 4: Linkage-method comparison

### SP4 Data model

This comparison is again carried out on structured synthetic data with known cluster labels ($n = 300$, $p = 500$;
cpp: $n = 1500$, $p = 50000$), so that different linkage rules can be judged by their ARI and runtime.

### SP4 Mathematical approach

The main question is how clustering changes when the dendrogram update rule changes. The demo compares several linkage
principles, including:

- **Single linkage**, based on nearest-neighbor inter-cluster distance,
- **Complete linkage**, based on farthest-neighbor inter-cluster distance,
- **Average linkage (UPGMA)**, based on average inter-cluster distance,
- **Weighted average linkage (WPGMA)**, based on weighted average inter-cluster distance,
- **Ward linkage**, based on minimum increase in within-cluster variance,
- and additional centroid/median-style variants, which don't have the reducibility property.

In abstract form, HAC iteratively merges two clusters $A$ and $B$ minimizing a linkage criterion

$$
d(A,B),
$$

where the definition of $d(A,B)$ depends on the chosen linkage rule. This subproblem studies how that choice affects
both accuracy and computational behavior on the same synthetic structure. The comparison here is run under the **exact**
Pearson-correlation distance, and it also exercises centroid- and median-style variants that lack the reducibility
property. (The C++ demo runs single linkage twice — once through the SLINK class directly and once through the
dispatcher; the Python binding has a single `agglomerative_cluster()` entry point, so single linkage appears once.)

---

## Subproblem 5: LSH-approximate linkage comparison

### SP5 Data model

This subproblem uses the same block-correlated latent-factor model as Subproblems 2–4, but on a larger, imbalanced
configuration ($n = 500$, $p = 1000$; cpp: $n = 50000$, $p = 100000$). Because the latent factors create dense,
well-separated blocks in correlation space, the geometry is favorable for approximate hashing.

### SP5 Mathematical approach

The comparison repeats the linkage study of Subproblem 4, but replaces the exact correlation distance with the pure
$O(1)$ SimHash-style LSH approximation

$$
d_{\mathrm{LSH}}(j,\ell)
\approx
1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
$$

Single (SLINK), complete, average (UPGMA), and weighted average (WPGMA) linkage are each run under this approximation.
The purpose is to show that, when the block geometry is well separated, the approximate distance still recovers the
correct structure at a substantially reduced per-distance cost.

---

## Subproblem 6: AR(1) Toeplitz block clustering

### SP6 Data model

Here the variables ($n = 300$, $p = 500$; cpp: $n = 3000$, $p = 50000$) are generated with a **block-diagonal AR(1)
Toeplitz** correlation structure. Each block $k$ is built by the stationary recursion
$\boldsymbol{x}_{j} = \rho_{k}\,\boldsymbol{x}_{j-1} + \sqrt{1 - \rho_{k}^{2}}\,\boldsymbol{\varepsilon}_{j}$ with
$\boldsymbol{\varepsilon}_{j}$ standard normal, which leaves every column at unit variance. The correlation therefore
decays *exactly* with the index separation inside a block, and vanishes between blocks:

$$
\mathrm{corr}(\boldsymbol{x}_{j}, \boldsymbol{x}_{j + m}) =
\begin{cases}
\rho_{k}^{m}, & b(j) = b(j+m) = k, \\
0, & \text{otherwise},
\end{cases}
$$

so a variable is strongly correlated with its immediate neighbors but nearly uncorrelated with distant variables in the
same block. Each block carries its own decay rate — $\rho_{k} = 0.95,\, 0.90,\, 0.85,\, 0.80,\, 0.75$ — so even the
slowest-decaying block is down to $0.95^{100} \approx 0.006$ after a hundred steps.

### SP6 Mathematical approach

This chaining geometry is used to contrast linkage rules under the **exact** Pearson-correlation distance. Single
linkage (SLINK) follows the correlation gradient and chains each block together, whereas complete and average linkage
tend to fragment the gradually decaying blocks. LSH is deliberately not used here, because the coarse angular
quantization destroys the continuous correlation gradient.

---

## Running the demo

```bash
python demo_mlm_hac_01.py
```

The demo prints its main diagnostics to the console.

---

## What to look for

When reading the output, focus on the following points:

- whether the recovered clusters match the known synthetic structure,
- how strongly ARI depends on the data geometry,
- whether the LSH approximation preserves the dominant block structure,
- how the linkage rule changes both runtime and clustering agreement — in particular, the AR(1) scenario should show
  single linkage at ARI $\approx 1$ while complete and average linkage shatter the chained blocks.

---

## Technical notes

- `agglomerative_cluster()` clusters the **rows** of its input and returns a SciPy-style $[N-1, 4]$ linkage matrix;
  `cut_tree()` recovers flat labels. Sample clustering passes `X` directly, variable clustering passes
  `np.asfortranarray(X.T)`. (The C++ engine clusters columns and uses Eigen transposed views instead.)
- All six subproblems are enabled in the main block; each is guarded by its own `if True:` toggle so an individual part
  can be switched off while iterating on another (mirroring the cpp `main()`).
- Subproblems 4 and 6 use the exact Pearson-correlation distance, while Subproblems 3 and 5 exercise the LSH filter
  and/or pure $O(1)$ LSH approximation.
- The cpp demo fills matrices with OpenMP-parallel per-thread RNGs (C++-exclusive); vectorized NumPy generation is used
  here, so numbers match the cpp output qualitatively, not bitwise.
- The C++ thread count (6) is mirrored via `set_num_threads(min(6, get_max_threads()))` from `trex_selector_neo.utils`.

---

**Last updated**: 2026-07-21
