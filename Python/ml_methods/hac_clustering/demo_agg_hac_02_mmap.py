# ==============================================================================
# demo_agg_hac_02_mmap.py
# ==============================================================================
#
# Demonstrates agglomerative hierarchical clustering on memory-mapped data.
# Mirrors R/ml_methods/hac_clustering/demo_agg_hac_02_mmap.R.
#
# Demo content:
#   - Converting an in-memory matrix to a memory-mapped file
#   - Clustering with different linkage methods on mmap-backed data
#   - Using different distance metrics
#   - Cutting the dendrogram to obtain cluster assignments
#
# Key API notes:
#   agglomerative_cluster() always takes np.ndarray for `data`.
#   X_mmap.to_numpy() provides a zero-copy view backed by the mmap file.
#   use_mmap=True controls the C++ backend's *internal* storage of
#   intermediate distance data — it does not change the data argument type.
#   The cluster assignments produced here are identical to demo_01 because
#   the same data is used; only the internal computation storage differs.
#
# ==============================================================================

import os
import tempfile

import numpy as np
from trex_selector.ml_methods.clustering import (
    agglomerative_cluster,
    cut_tree,
    DistanceMetric,
    LinkageMethod,
)
from trex_selector.utils import numpy_to_memmap

# ==============================================================================
# Global Parameters and Data Generation
# ==============================================================================

rng = np.random.default_rng(42)

n = 50
p = 10
X = rng.standard_normal((n, p))

num_clusters = 5

# Create temporary mmap backing file
fd, mmap_path = tempfile.mkstemp(suffix=".bin")
os.close(fd)

print("\n" + "=" * 70)
print("Agglomerative Hierarchical Clustering (memory-mapped)")
print("=" * 70 + "\n")
print(f"Data: n = {n} samples, p = {p} features (mmap-backed)\n")

try:
    # Write X to the mmap file and obtain a zero-copy numpy view.
    # X_view is a Fortran-contiguous float64 ndarray backed by the binary file;
    # agglomerative_cluster reads from it without pulling the entire matrix into
    # a separate Python buffer.
    X_mmap = numpy_to_memmap(mmap_path, X)
    X_view = X_mmap.to_numpy()

    # ==========================================================================
    # Part A: Ward linkage with Euclidean distance
    # ==========================================================================

    print("Part A: Ward linkage with Euclidean distance (mmap)")
    print("-" * 70 + "\n")

    linkage_ward_euclidean = agglomerative_cluster(
        X_view,
        method=LinkageMethod.Ward,
        metric=DistanceMetric.Euclidean,
        use_mmap=True,
    )

    print(f"Linkage matrix shape: {linkage_ward_euclidean.shape}")
    print("First 5 rows of linkage matrix:")
    print(np.round(linkage_ward_euclidean[:5], 4))

    clusters_ward = cut_tree(
        linkage_ward_euclidean,
        num_orig_objs=n,
        num_clusters=num_clusters,
    )
    print(f"\nCluster assignments (k = {num_clusters}):")
    print(f"Cluster distribution: {np.bincount(clusters_ward).tolist()}\n")

    # ==========================================================================
    # Part B: Average linkage with Euclidean distance
    # ==========================================================================

    print("Part B: Average linkage with Euclidean distance (mmap)")
    print("-" * 70 + "\n")

    linkage_avg_euclidean = agglomerative_cluster(
        X_view,
        method=LinkageMethod.Average,
        metric=DistanceMetric.Euclidean,
        use_mmap=True,
    )

    clusters_avg = cut_tree(
        linkage_avg_euclidean,
        num_orig_objs=n,
        num_clusters=num_clusters,
    )
    print(f"Cluster assignments (k = {num_clusters}):")
    print(f"Cluster distribution: {np.bincount(clusters_avg).tolist()}\n")

    # ==========================================================================
    # Part C: Complete linkage with Euclidean distance
    # ==========================================================================

    print("Part C: Complete linkage with Euclidean distance (mmap)")
    print("-" * 70 + "\n")

    linkage_comp_euclidean = agglomerative_cluster(
        X_view,
        method=LinkageMethod.Complete,
        metric=DistanceMetric.Euclidean,
        use_mmap=True,
    )

    clusters_comp = cut_tree(
        linkage_comp_euclidean,
        num_orig_objs=n,
        num_clusters=num_clusters,
    )
    print(f"Cluster assignments (k = {num_clusters}):")
    print(f"Cluster distribution: {np.bincount(clusters_comp).tolist()}\n")

    # ==========================================================================
    # Part D: Ward linkage with Correlation distance
    # ==========================================================================

    print("Part D: Ward linkage with Correlation distance (mmap)")
    print("-" * 70 + "\n")

    linkage_ward_corr = agglomerative_cluster(
        X_view,
        method=LinkageMethod.Ward,
        metric=DistanceMetric.Correlation,
        use_mmap=True,
    )

    clusters_ward_corr = cut_tree(
        linkage_ward_corr,
        num_orig_objs=n,
        num_clusters=num_clusters,
    )
    print(f"Cluster assignments (k = {num_clusters}):")
    print(f"Cluster distribution: {np.bincount(clusters_ward_corr).tolist()}\n")

finally:
    try:
        os.unlink(mmap_path)
    except OSError:
        pass

print("=" * 70)
print("demo_agg_hac_02_mmap.py complete.")
