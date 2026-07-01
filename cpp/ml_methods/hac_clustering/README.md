# HAC Clustering: Demonstration Suite

## Overview

This folder contains C++ examples for **hierarchical agglomerative clustering (HAC)** in the `ml_methods` module.

Within the TRexSelector project, HAC is used by the **GVS (Grouped Variable Selection)** workflow to construct variable groups from a dendrogram cut. The examples in this folder show how clustering behaves in smaller in-memory settings and how the same ideas can be scaled to much larger problems through memory-mapped I/O.

The main goals of this folder are:

1. to demonstrate how HAC is applied to synthetic sample and variable clustering tasks,
2. to compare clustering behavior under different distance and linkage choices,
3. to illustrate how large clustering problems can be handled with memory-mapped storage,
4. to provide a bridge between exploratory clustering demos and downstream grouped-variable workflows.

---

## What this folder covers

Hierarchical agglomerative clustering builds a tree of successive merges and then recovers cluster assignments by cutting that tree at a chosen level.

In this folder, the demos focus on:

- **in-memory HAC experiments**, for understanding clustering structure and linkage behavior,
- **correlation-based variable clustering**, which is relevant for grouped variable selection,
- **memory-mapped large-scale HAC**, for problems that would otherwise exceed available RAM.

---

## Start here

If you are new to this folder, begin with:

1. **`demo_mlm_hac_01/`**  
   An in-memory HAC demo suite covering several clustering scenarios and linkage comparisons.

2. **`demo_mlm_hac_02_mmap/`**  
   A memory-mapped large-scale demo that reproduces the LSH linkage comparison idea at biobank-representative scale.

---

## Demo overview

| Folder                   | Purpose |
|--------------------------|---------|
| `demo_mlm_hac_01/`       | In-memory HAC demo suite for sample clustering, variable clustering, and linkage-method comparisons |
| `demo_mlm_hac_02_mmap/`  | Memory-mapped large-scale HAC demo for LSH-based linkage comparison without full in-RAM allocation |

---

## Running

```bash
./build/debug/bin/ml_methods/hac_clustering/demo_mlm_hac_01/demo_mlm_hac_01
./build/debug/bin/ml_methods/hac_clustering/demo_mlm_hac_02_mmap/demo_mlm_hac_02_mmap
```

---

## Folder contents

```txt
hac_clustering/
  ├── README.md
  ├── demo_mlm_hac_01/
  │   ├── demo_mlm_hac_01.cpp
  │   └── README.md
  └── demo_mlm_hac_02_mmap/
      ├── demo_mlm_hac_02_mmap.cpp
      └── README.md
```

---

## Notes for new users

- `demo_mlm_hac_01` is best for understanding the clustering logic in standard in-memory workflows.
- `demo_mlm_hac_02_mmap` is intended for very large synthetic settings where the matrix should not be allocated fully in RAM.
- Both demos print their main results to the console.
- Validation-oriented cross-checks are linked from the demo-specific README files where relevant.

---

**Last updated**: 2026-07-01
