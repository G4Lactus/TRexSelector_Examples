# HAC Clustering: Demonstration Suite

## Overview

This folder contains Python examples for **hierarchical agglomerative clustering (HAC)** in the
`trex_selector_neo.ml_methods.clustering` module.

Within the TRexSelector project, HAC is used by the **GVS (Grouped Variable Selection)** workflow to construct variable
groups from a dendrogram cut. The examples in this folder show how clustering behaves in in-memory settings and how the
same ideas extend to memory-mapped, out-of-core workflows.

The main goals of this folder are:

1. to demonstrate how HAC is applied to synthetic sample and variable clustering tasks,
2. to compare clustering behavior under different distance and linkage choices,
3. to illustrate how large clustering problems can be handled with memory-mapped storage,
4. to provide a bridge between exploratory clustering demos and downstream grouped-variable workflows.

The C++ counterparts live in [cpp/ml_methods/hac_clustering/](../../../cpp/ml_methods/hac_clustering/); equivalent R
walkthroughs are available under `R/ml_methods/hac_clustering/`.

---

## What this folder covers

Hierarchical agglomerative clustering builds a tree of successive merges and then recovers cluster assignments by
cutting that tree at a chosen level.

In this folder, the demos focus on:

- **in-memory HAC experiments**, for understanding clustering structure and linkage behavior,
- **correlation-based variable clustering**, which is relevant for grouped variable selection,
- **memory-mapped HAC**, for problems that would otherwise exceed available RAM.

The public entry points are `agglomerative_cluster(X, method, metric, use_mmap)` — which clusters the **rows** of `X`
and returns a SciPy-style linkage matrix — and `cut_tree(linkage, num_orig_objs, num_clusters)`, together with the
`LinkageMethod` and `DistanceMetric` enums.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_hac_01/`**  
   An in-memory HAC demo suite covering several clustering scenarios and linkage comparisons.

2. **`demo_mlm_hac_02_mmap/`**  
   A memory-mapped demo that reproduces the C++ biobank-scale LSH linkage comparison at an illustrative scale.

---

## Demo overview

| Folder | Purpose |
| -------------------------- | --------- |
| [demo_mlm_hac_01/](demo_mlm_hac_01/README.md) | In-memory HAC demo suite for sample clustering, variable clustering, LSH benchmarks, and linkage-method comparisons |
| [demo_mlm_hac_02_mmap/](demo_mlm_hac_02_mmap/README.md) | Memory-mapped HAC demo for LSH-based linkage comparison without full in-RAM allocation |

---

## Running

```bash
python demo_mlm_hac_01/demo_mlm_hac_01.py
python demo_mlm_hac_02_mmap/demo_mlm_hac_02_mmap.py
```

---

## Folder contents

```txt
hac_clustering/
  ├── README.md
  ├── demo_mlm_hac_01/
  │   ├── demo_mlm_hac_01.py
  │   └── README.md
  └── demo_mlm_hac_02_mmap/
      ├── demo_mlm_hac_02_mmap.py
      └── README.md
```

---

## Notes for new users

- `demo_mlm_hac_01` is best for understanding the clustering logic in standard in-memory workflows.
- `demo_mlm_hac_02_mmap` shows the out-of-core workflow; the C++ counterpart runs it at ~298 GB scale, the Python
  mirror at ~61 MB.
- Both demos run at sizes scaled down from their C++ counterparts (noted in the script headers) so that they finish in
  well under a second.
- Both demos print their main results to the console.
- Cross-language correctness checks against R references live in the TRexSelector library test suite
  (`cpp/tests/validation/`), run on demand rather than shipped with these examples.

---

**Last updated**: 2026-07-21
