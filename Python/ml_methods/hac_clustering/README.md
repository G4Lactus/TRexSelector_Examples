# HAC Clustering Demos (Python)

## Purpose

Demonstrate hierarchical agglomerative clustering (HAC) via the
`trex_selector_neo.ml_methods.clustering` module — building a linkage matrix
with different linkage methods and distance metrics, and cutting the
dendrogram into cluster assignments.

APIs used:

- `agglomerative_cluster(X, method=..., metric=..., use_mmap=...)` — returns
  the linkage matrix
- `cut_tree(linkage, num_orig_objs=..., num_clusters=...)` — returns cluster
  assignments
- `LinkageMethod` — `Ward`, `Average`, `Complete`
- `DistanceMetric` — `Euclidean`, `Correlation`
- `numpy_to_memmap` (from `trex_selector_neo.utils`, mmap variant only)

---

## Demos

| Demo | Description |
|---|---|
| [demo_agg_hac_01_in_memory.py](demo_agg_hac_01_in_memory.py) | HAC on an in-memory 50 x 10 Gaussian matrix (`use_mmap=False`). Parts A–D: Ward/Euclidean, Average/Euclidean, Complete/Euclidean, Average/Correlation (Ward is Euclidean-only); each linkage is cut into k = 5 clusters via `cut_tree`. |
| [demo_agg_hac_02_mmap.py](demo_agg_hac_02_mmap.py) | Same four linkage/metric combinations, but the data is first written to a temporary binary file via `numpy_to_memmap()` and passed as the zero-copy `to_numpy()` view, with `use_mmap=True` so the C++ backend also stores intermediate distance data on disk. |

Key API note (from the mmap demo): `agglomerative_cluster()` always takes an
`np.ndarray` for `data`; `use_mmap=True` only controls the backend's
*internal* storage of intermediate distance data. Since both demos use the
same data and seed, the cluster assignments are identical — only the internal
computation storage differs.

---

## Running the Demos

```bash
python demo_agg_hac_01_in_memory.py
python demo_agg_hac_02_mmap.py
```

Output is printed to the console (linkage-matrix excerpts and cluster-size
distributions). The mmap demo cleans up its temporary binary file on exit.

---

## Counterparts

- C++: [cpp/ml_methods/hac_clustering/demo_mlm_hac_01](../../../cpp/ml_methods/hac_clustering/)
  and [demo_mlm_hac_02_mmap](../../../cpp/ml_methods/hac_clustering/)
- R: `R/ml_methods/hac_clustering/demo_agg_hac_01_in_memory.R` and
  `demo_agg_hac_02_mmap.R` (these scripts are direct mirrors)

---

**Last updated**: 2026-07-06
