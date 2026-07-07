# ==============================================================================
# demo_agg_hac_01_in_memory.py
# ==============================================================================
#
# Demonstrates agglomerative hierarchical clustering on in-memory data.
# Mirrors R/ml_methods/hac_clustering/demo_agg_hac_01_in_memory.R.
#
# Demo content:
#   - Clustering with different linkage methods (Ward, Average, Complete)
#   - Using different distance metrics (Euclidean, Correlation)
#   - Cutting the dendrogram to obtain cluster assignments
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods.clustering import (
    agglomerative_cluster,
    cut_tree,
    DistanceMetric,
    LinkageMethod,
)

# ==============================================================================
# Global Parameters and Data Generation
# ==============================================================================

rng = np.random.default_rng(42)

n = 50
p = 10
X = rng.standard_normal((n, p))

num_clusters = 5

print("\n" + "=" * 70)
print("Agglomerative Hierarchical Clustering (in-memory)")
print("=" * 70 + "\n")
print(f"Data: n = {n} samples, p = {p} features\n")


# ==============================================================================
# Part A: Ward linkage with Euclidean distance
# ==============================================================================

print("Part A: Ward linkage with Euclidean distance")
print("-" * 70 + "\n")

linkage_ward_euclidean = agglomerative_cluster(
    X,
    method=LinkageMethod.Ward,
    metric=DistanceMetric.Euclidean,
    use_mmap=False,
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


# ==============================================================================
# Part B: Average linkage with Euclidean distance
# ==============================================================================

print("Part B: Average linkage with Euclidean distance")
print("-" * 70 + "\n")

linkage_avg_euclidean = agglomerative_cluster(
    X,
    method=LinkageMethod.Average,
    metric=DistanceMetric.Euclidean,
    use_mmap=False,
)

clusters_avg = cut_tree(
    linkage_avg_euclidean,
    num_orig_objs=n,
    num_clusters=num_clusters,
)
print(f"Cluster assignments (k = {num_clusters}):")
print(f"Cluster distribution: {np.bincount(clusters_avg).tolist()}\n")


# ==============================================================================
# Part C: Complete linkage with Euclidean distance
# ==============================================================================

print("Part C: Complete linkage with Euclidean distance")
print("-" * 70 + "\n")

linkage_comp_euclidean = agglomerative_cluster(
    X,
    method=LinkageMethod.Complete,
    metric=DistanceMetric.Euclidean,
    use_mmap=False,
)

clusters_comp = cut_tree(
    linkage_comp_euclidean,
    num_orig_objs=n,
    num_clusters=num_clusters,
)
print(f"Cluster assignments (k = {num_clusters}):")
print(f"Cluster distribution: {np.bincount(clusters_comp).tolist()}\n")


# ==============================================================================
# Part D: Average linkage with Correlation distance
# ==============================================================================
# Ward (like Median/Centroid) is only valid in Euclidean geometry — the library
# rejects Ward with non-Euclidean metrics. Graph-based linkages such as Average
# or Complete work with any metric, so we pair Correlation with Average here.

print("Part D: Average linkage with Correlation distance")
print("-" * 70 + "\n")

linkage_avg_corr = agglomerative_cluster(
    X,
    method=LinkageMethod.Average,
    metric=DistanceMetric.Correlation,
    use_mmap=False,
)

clusters_avg_corr = cut_tree(
    linkage_avg_corr,
    num_orig_objs=n,
    num_clusters=num_clusters,
)
print(f"Cluster assignments (k = {num_clusters}):")
print(f"Cluster distribution: {np.bincount(clusters_avg_corr).tolist()}\n")


print("=" * 70)
print("demo_agg_hac_01_in_memory.py complete.")
