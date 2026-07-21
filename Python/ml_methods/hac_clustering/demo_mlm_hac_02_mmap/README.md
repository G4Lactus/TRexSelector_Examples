# Demo: Large-Scale HAC with Memory-Mapped I/O

## Purpose

This demo illustrates how hierarchical agglomerative clustering (HAC) can be applied using memory-mapped storage.

It extends the block-correlated variable-clustering setting from the in-memory HAC examples to an out-of-core workflow
where the data matrix lives in a file rather than in RAM. The demo combines synthetic data generation directly into the
mapped file, in-place preprocessing, and clustering under an LSH-based correlation approximation.

The workflow is organized into three phases:

1. **Generation**: write a correlated matrix directly into a memory-mapped file.
2. **Standardization**: mean-center and L2-normalize the variables in place.
3. **Clustering**: benchmark selected HAC linkage strategies under an LSH-based correlation distance approximation.

It mirrors the corresponding C++ example in
[cpp/ml_methods/hac_clustering/demo_mlm_hac_02_mmap/](../../../../cpp/ml_methods/hac_clustering/demo_mlm_hac_02_mmap/),
which runs at biobank-representative scale ($n = 100000$, $p = 400000$, ~298 GB on disk). The Python mirror keeps the
same three-phase structure, block layout, and test sequence at an illustrative scale (~61 MB on disk).

---

## Statistical model and notation

Let

- $n$ the number of samples,
- $p$ the number of variables,
- $K$ the number of true variable blocks,
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ denote the data matrix,
- $b(j) \in \lbrace1,\dots,K\rbrace$ the true block label of variable $j$.

In this demo, the dimensions are

$$
n = 2000, \qquad p = 4000, \qquad K = 10
$$

(cpp: $n = 100000$, $p = 400000$, $K = 10$; the ten block sizes are the cpp block sizes scaled by $1/100$).

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

The full matrix is not assembled in RAM. Instead, the entries of $\boldsymbol{X}$ are written variable by variable into
a memory-mapped file, while only the latent signal matrix

$$
\boldsymbol{L} \in \mathbb{R}^{n \times K}
$$

is held in memory.

This changes the storage pattern from a RAM-bound dense allocation to an out-of-core workflow. Because the Python
binding clusters the **rows** of its input, the matrix is stored **transposed** relative to the C++ demo: the file
holds a $p \times n$ layout with one variable per row, so the zero-copy `to_numpy()` view can be passed to the
clustering call directly.

---

## Phase 2: In-place standardization

Before correlation-based clustering, each variable is standardized in place by

$$
\widetilde{\boldsymbol{x}}_j =
\frac{\boldsymbol{x}_j - \bar{x}_j \mathbf{1}}
{\lVert \boldsymbol{x}_j - \bar{x}_j \mathbf{1} \rVert_2}.
$$

This preprocessing step is required for correlation-based distance calculations and is performed directly on the
memory-mapped matrix (one variable at a time, sequential-streaming-friendly) rather than on a copied in-memory version.

---

## Phase 3: Clustering

After standardization, the demo clusters variables under the **LSH correlation approximation**,

$$
d_{\mathrm{LSH}}(j,\ell)
\approx
1 - \left| \mathrm{corr}(\boldsymbol{x}_j, \boldsymbol{x}_\ell) \right|,
$$

a pure $O(1)$ SimHash-based correlation approximation, used to avoid repeatedly evaluating full high-dimensional exact
correlations at scale. The `use_mmap=True` flag additionally keeps the backend's intermediate distance data on disk.

### Subproblem 1: Single linkage

The first clustering task applies **single linkage** under the LSH approximation. Algorithmically, this is dispatched
to the single-linkage HAC path and evaluated by cutting the dendrogram into $K$ groups and comparing the result against
the known block labels via ARI.

### Subproblem 2: Complete linkage and WPGMA feasibility

The second subproblem is a **feasibility question** rather than a benchmark run. Complete linkage and WPGMA are skipped
because they require the $O(p^2)$ Lance-Williams matrix machinery and are therefore intractable at the C++ demo's scale
($p = 400000$, ~1.2 TB of NVMe I/O). Although they would be tractable at the reduced Python scale, the skip is kept to
mirror the C++ console output and its reasoning.

### Subproblem 3: Average linkage

The third subproblem applies **average linkage (UPGMA)** under the same LSH-based distance approximation. The backend
routes this through a projected low-dimensional internal representation rather than performing the full clustering
directly on a quadratic out-of-core structure.

### Subproblem 4: Ward linkage

The fourth subproblem applies **Ward’s minimum-variance linkage** under the same setup. As in the average-linkage case,
the backend uses an internal projected engine so that the clustering remains computationally manageable at high
dimension.

---

## Running the demo

```bash
python demo_mlm_hac_02_mmap.py
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

You can compare this out-of-core behavior with the in-memory HAC examples in
[demo_mlm_hac_01](../demo_mlm_hac_01/README.md).

---

## Technical notes

- The C++ demo injects the file path at compile time via `TREX_DEMO_MMAP_PATH` and skips Phases 1–2 when a file of the
  expected size already exists; the Python mirror creates a fresh temporary file with `tempfile.mkstemp()` on every run
  and always regenerates.
- The C++ demo generates and standardizes with OpenMP parallel loops (C++-exclusive); the Python mirror streams one
  variable at a time through the zero-copy view, so numbers match qualitatively, not bitwise.
- Cleanup is handled in a `try/finally` block: the mapped objects are deleted (unmapped) first, then the backing file
  is removed.

---

**Last updated**: 2026-07-21
