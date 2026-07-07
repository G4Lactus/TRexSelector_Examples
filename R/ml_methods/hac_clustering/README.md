# HAC Clustering — R Demos

## Purpose

Demonstrates hierarchical agglomerative clustering (HAC) through the
**TRexSelectorNeo** function `agglomerative_cluster()`, on both in-memory and
memory-mapped data. Every scenario reports the Adjusted Rand Index (ARI)
against known ground-truth cluster labels.

---

## Demos

| File | Description |
|---|---|
| [demo_agg_hac_01_in_memory.R](demo_agg_hac_01_in_memory.R) | HAC on in-memory data, five scenarios with ARI evaluation (mirrors C++ `demo_mlm_hac_01.cpp`) |
| [demo_agg_hac_02_mmap.R](demo_agg_hac_02_mmap.R) | Same clustering routed through mmap-backed storage (`use_mmap = TRUE`); confirms ARI is identical via the mmap path |
| [hac_helpers.R](hac_helpers.R) | Shared helpers sourced by both demos — not run standalone |

### Demo 01 scenarios

1. **Sample clustering (Euclidean)** — structured row data with k=3 Gaussian
   cluster centers; Euclidean-geometry linkage set (Ward, Average, Complete,
   Single, WPGMA, Median, Centroid).
2. **Variable clustering (Correlation)** — latent-factor block-correlated
   column data with k=5 blocks; general linkage set (Average, Complete,
   Single, WPGMA — Ward/Median/Centroid excluded as they require Euclidean
   geometry).
3. **Sample clustering (Manhattan)** — same row data as Scenario 1; general
   linkage set.
4. **Variable clustering (Correlation_LSH_Filter)** — exact correlation when
   the Hamming distance is <= 14, otherwise 1.0.
5. **Variable clustering (Correlation_LSH_Approx)** — O(1) SimHash
   approximation `1 - |cos(pi * hamming / 64)|`; Ward included via a 64-D
   SimHash-projected Euclidean space (the one non-Euclidean metric for which
   Ward is a defined, approximate operation).

Demo 02 repeats Scenarios 1 and 2 with the design (or its transpose, for
variable clustering) written to an mmap file first.

### hac_helpers.R contents

- `compute_ari()` — Adjusted Rand Index, dependency-free R port of the inline
  C++ implementation in `demo_mlm_hac_01.cpp`.
- `generate_sample_clusters()` — structured row data with distinct Gaussian
  centers.
- `generate_variable_blocks()` — latent-factor block-correlated column data.
- `standardize_for_correlation()` — mean-center + L2-normalize each column.

### APIs used

`agglomerative_cluster()`, the `LinkageMethod$...` enum (Ward, Average,
Complete, Single, WPGMA, Median, Centroid), and — in demo 02 —
`convert_to_memory_mapped()` for the mmap path.

---

## Running

```bash
Rscript R/ml_methods/hac_clustering/demo_agg_hac_01_in_memory.R
Rscript R/ml_methods/hac_clustering/demo_agg_hac_02_mmap.R
```

Console output only.

---

## Counterparts

- C++: [cpp/ml_methods/hac_clustering/](../../../cpp/ml_methods/hac_clustering/)
  — `demo_mlm_hac_01` and `demo_mlm_hac_02_mmap`.
- Cross-validation programs against R reference runs live in
  [cpp/ml_methods/validation/hac_clustering/](../../../cpp/ml_methods/validation/).

---

**Last updated**: 2026-07-06
