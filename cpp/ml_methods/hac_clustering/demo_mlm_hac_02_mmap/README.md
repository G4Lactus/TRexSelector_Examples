# Demo: Large-Scale HAC with Memory-Mapped I/O

## Purpose

This demo illustrates how hierarchical agglomerative clustering (HAC) can be applied using memory-mapped storage.

It extends the block-correlated variable-clustering setting from the in-memory HAC examples to a setting where the full
matrix is too large to allocate in RAM. The demo combines large-scale synthetic data generation, in-place preprocessing,
and clustering under an LSH-based correlation approximation.

The workflow is organized into three phases:

1. **Generation**: write a large correlated matrix directly into a memory-mapped file.
2. **Standardization**: mean-center and L2-normalize columns in place.
3. **Clustering**: benchmark selected HAC linkage strategies under an LSH-based correlation distance approximation.

---

## Statistical model and notation

Let

- $n$ the number of samples,
- $p$ the number of variables,
- $K$ the number of true variable blocks,
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $b(j) \in \{1,\dots,K\}$ the true block label of variable $j$.

In this demo, the dimensions are

$$
n = 100000, \qquad p = 400000, \qquad K = 10.
$$

The quality of the recovered clustering is assessed with the Adjusted Rand Index (ARI), where $\mathrm{ARI} = 1$
indicates perfect recovery of the true block structure.

---

## Data model

The synthetic variables are generated from a blockwise latent-factor model. For each variable $j$ in block $b(j)$, the
entries of $\boldsymbol{X}$ are generated according to

$$
X_{ij} =
\mu_{b(j)}
+
\sigma_{b(j)}
\left[
\alpha L_{i,b(j)} + \gamma \varepsilon_{ij}
\right],
$$

where

- $L_{i,b(j)} \sim \mathcal{N}(0,1)$ is the latent signal shared within block $b(j)$,
- $\varepsilon_{ij} \sim \mathcal{N}(0,1)$ is independent noise,
- $\mu_{b(j)}$ is a block-specific mean,
- $\sigma_{b(j)}$ is a block-specific scale.

In the implementation, the within-block signal and noise weights are fixed to

$$
\alpha = 0.85,
\qquad
\gamma = 0.15,
$$

which produces strongly correlated variable blocks while preserving some within-block heterogeneity.

The true labels are determined deterministically from the block sizes, so the clustering target is known exactly before
the HAC stage begins.

---

## Phase 1: Memory-mapped generation

The full matrix is not assembled in RAM. Instead, the entries of $\boldsymbol{X}$ are written column by column into a
memory-mapped file, while only the latent signal matrix

$$
\boldsymbol{L} \in \mathbb{R}^{n \times K}
$$

is held in memory.

This changes the storage pattern from a RAM-bound dense allocation to an out-of-core workflow. In the current
implementation, the memory-mapped file has a disk footprint of about 298 GB, whereas the latent $n \times K$ signal
matrix occupies only a small RAM budget relative to the full design matrix.

---

## Phase 2: In-place standardization

Before correlation-based clustering, each variable is standardized in place by

$$
\widetilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_2}.
$$

This preprocessing step is required for correlation-based distance calculations and is performed directly on the
memory-mapped matrix rather than on a copied in-memory version.

---

## Phase 3: Clustering

After standardization, the demo clusters variables under the **LSH correlation approximation**,

$$
d_{\mathrm{LSH}}(j,\ell)
\approx
1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|.
$$

The source describes this as a pure $O(1)$ SimHash-based correlation approximation, used to avoid repeatedly evaluating
full high-dimensional exact correlations at this scale.

### Subproblem 1: Single linkage

The first clustering task applies **single linkage** under the LSH approximation. Algorithmically, this is dispatched to
the single-linkage HAC path and evaluated by cutting the dendrogram into $K$ groups and comparing the result against the
known block labels via ARI.

### Subproblem 2: Complete linkage and WPGMA feasibility

The second subproblem is a **feasibility question** rather than a full benchmark run. In the current implementation,
complete linkage and WPGMA are skipped because they require the $O(p^2)$ Lance-Williams matrix machinery and are
therefore treated as computationally intractable at $p = 400000$ in this setting.

### Subproblem 3: Average linkage

The third subproblem applies **average linkage (UPGMA)** under the same LSH-based distance approximation. The
implementation routes this through a projected low-dimensional internal representation rather than performing the full
clustering directly on a quadratic out-of-core structure.

### Subproblem 4: Ward linkage

The fourth subproblem applies **Ward’s minimum-variance linkage** under the same large-scale setup. As in the
average-linkage case, the implementation uses an internal projected engine so that the clustering remains
computationally manageable at this dimension.

---

## Running the demo

```bash
./build/debug/bin/ml_methods/hac_clustering/demo_mlm_hac_02_mmap/demo_mlm_hac_02_mmap
```

The demo prints its diagnostics to the console, including file size, latent-memory footprint, phase timings, and ARI
values for the executed clustering tasks.

---

## What to look for

When reading the output, focus on the following points:

- whether ARI remains close to 1 despite the use of approximate correlation distances,
- whether the block structure is recovered consistently across the executed linkage strategies,
- how the three workflow phases differ in runtime,
- which linkage strategies are executed and which are skipped for computational reasons,
- whether the temporary memory-mapped file is cleaned up successfully at the end.

You can compare the large-scale behavior with the in-memory HAC examples in
[demo_mlm_hac_01](../demo_mlm_hac_01/README.md).

The same memory-mapped clustering technique is demonstrated at a smaller, illustrative scale in
[Python](../../../../Python/ml_methods/hac_clustering/demo_agg_hac_02_mmap.py) and
[R](../../../../R/ml_methods/hac_clustering/demo_agg_hac_02_mmap.R).

---

## Technical notes

- The memory-mapped file path is injected at compile time via `TREX_DEMO_MMAP_PATH`, with a relative fallback if not
  provided.
- If a file with the expected size already exists, the generation and standardization phases are skipped and the demo
  reuses the existing mapped matrix for clustering. Because the file is pre-allocated to its full size before any data
  is written, a run that is interrupted mid-generation can leave behind a same-sized but incomplete file that this
  check cannot distinguish from a valid one.
- Cleanup is handled automatically on normal exit by unmapping first and then removing the backing file.

---

**Last updated**: 2026-07-09
