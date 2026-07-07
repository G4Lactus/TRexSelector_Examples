# ==============================================================================
# demo_agg_hac_01_in_memory.R
# ==============================================================================
#
# The file demonstrates agglomerative hierarchical clustering on in-memory data
# with Adjusted Rand Index (ARI) performance evaluation.
#
# Scenario 1: Sample Clustering (Euclidean)
#   Structured row data with k=3 distinct Gaussian cluster centers.
#   Euclidean-geometry linkage set: Ward, Average, Complete, Single, WPGMA,
#   Median, Centroid.
#   Mirrors C++ demo_mlm_hac_01.cpp — Demo 1.
#
# Scenario 2: Variable Clustering (Correlation)
#   Latent-factor block-correlated column data with k=5 blocks.
#   General linkage set: Average, Complete, Single, WPGMA.
#   (Ward/Median/Centroid require Euclidean geometry and are excluded.)
#   Mirrors C++ demo_mlm_hac_01.cpp — Demo 2.
#
# Scenario 3: Sample Clustering (Manhattan)
#   Same row-data setup as Scenario 1; general linkage set.
#
# Scenario 4: Variable Clustering (Correlation_LSH_Filter)
#   Same block-correlated setup as Scenario 2; general linkage set.
#   Exact correlation when Hamming distance <= 14, otherwise 1.0.
#
# Scenario 5: Variable Clustering (Correlation_LSH_Approx)
#   Same block-correlated setup as Scenario 2.
#   Ward included: routed via 64D SimHash-projected Euclidean space.
#   O(1) approximation: 1 - |cos(pi * hamming / 64)|.
#
# Each scenario reports the Adjusted Rand Index (ARI) against known ground truth.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

source(file.path(this_dir_, "hac_helpers.R"))

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Agglomerative Hierarchical Clustering (in-memory)\n")
cat(strrep("=", 70), "\n\n")

# Linkage methods: Euclidean-geometry set (Ward / Median / Centroid require
# Euclidean distances — centroid-based Lance-Williams update is only valid in
# Euclidean space). Used for Scenario 1 only.
linkage_methods_eucl <- list(
  list(name = "Ward",             method = LinkageMethod$Ward),
  list(name = "Average (UPGMA)",  method = LinkageMethod$Average),
  list(name = "Complete",         method = LinkageMethod$Complete),
  list(name = "Single (SLINK)",   method = LinkageMethod$Single),
  list(name = "WPGMA",            method = LinkageMethod$WPGMA),
  list(name = "Median (WPGMC)",   method = LinkageMethod$Median),
  list(name = "Centroid (UPGMC)", method = LinkageMethod$Centroid)
)

# Linkage methods: general set — purely distance-based, metric-agnostic.
# Used for Scenarios 2, 3, 4.
linkage_methods_gen <- list(
  list(name = "Average (UPGMA)", method = LinkageMethod$Average),
  list(name = "Complete",        method = LinkageMethod$Complete),
  list(name = "Single (SLINK)",  method = LinkageMethod$Single),
  list(name = "WPGMA",           method = LinkageMethod$WPGMA)
)

# Linkage methods: LSH_Approx set — Ward is additionally permitted because the
# dispatcher routes it to ProjectedGeometricUpdatePolicy (Ward in 64D SimHash-
# projected Euclidean space), making it a defined approximate operation.
# Used for Scenario 5 only.
linkage_methods_lsh_approx <- list(
  list(name = "Ward (SimHash)",  method = LinkageMethod$Ward),
  list(name = "Average (UPGMA)", method = LinkageMethod$Average),
  list(name = "Complete",        method = LinkageMethod$Complete),
  list(name = "Single (SLINK)",  method = LinkageMethod$Single),
  list(name = "WPGMA",           method = LinkageMethod$WPGMA)
)


# ==============================================================================
# Scenario 1: Sample Clustering (Euclidean)
# ==============================================================================
# n=300 samples in k=3 equal clusters; centers at 0, 5, 20 (well separated).
# Mirrors C++ generate_hierarchical_row_data(): N=1500, P=50000, K=3.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("Scenario 1: Sample Clustering (Euclidean)\n")
cat(strrep("=", 70), "\n\n")

n1        <- 300L
p1        <- 5000L
k1        <- 3L
centers1  <- c(0.0, 5.0, 20.0)
std_devs1 <- c(1.0, 1.0, 1.0)

cat(sprintf("Generating data: n = %d, p = %d, k = %d ...\n", n1, p1, k1))
dat1 <- generate_sample_clusters(n1, p1, k1, centers1, std_devs1, seed = 42L)
cat("Done.\n\n")

cat(sprintf("%-22s  ARI\n", "Linkage Method"))
cat(strrep("-", 37), "\n")

for (lm in linkage_methods_eucl) {

  linkage_mat <- agglomerative_cluster(
    data     = dat1$X,
    method   = lm$method,
    metric   = DistanceMetric$Euclidean,
    use_mmap = FALSE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = n1, num_clusters = k1)
  ari         <- compute_ari(dat1$true_labels, pred_labels)
  cat(sprintf("%-22s  %.4f\n", lm$name, ari))
}


# ==============================================================================
# Scenario 2: Variable Clustering (Correlation)
# ==============================================================================
# n=1000 observations; p=5000 variables in k=5 latent-factor blocks.
# Standardize columns (mean=0, L2=1) before clustering.
# Pass t(X) so rows = variables; use metric=Correlation.
# Mirrors C++ generate_correlated_variables(): N=10000, P=50000, K=5.
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Scenario 2: Variable Clustering (Correlation)\n")
cat(strrep("=", 70), "\n\n")

n2           <- 1000L
p2           <- 5000L
k2           <- 5L
block_sizes2 <- c(2000L, 1000L, 800L, 700L, 500L)

cat(sprintf("Generating data: n = %d, p = %d, k = %d ...\n", n2, p2, k2))
dat2 <- generate_variable_blocks(n2, p2, k2, block_sizes2,
                                 signal_strength = 0.85,
                                 noise_strength  = 0.15,
                                 seed = 42L)
cat("Standardizing columns for correlation ...\n")
X2_std <- standardize_for_correlation(dat2$X)
X2_t   <- t(X2_std)   # p2 x n2 — rows are the variables to cluster
cat("Done.\n\n")

cat(sprintf("%-22s  ARI\n", "Linkage Method"))
cat(strrep("-", 37), "\n")

for (lm in linkage_methods_gen) {

  linkage_mat <- agglomerative_cluster(
    data     = X2_t,
    method   = lm$method,
    metric   = DistanceMetric$Correlation,
    use_mmap = FALSE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = p2, num_clusters = k2)
  ari         <- compute_ari(dat2$true_labels, pred_labels)
  cat(sprintf("%-22s  %.4f\n", lm$name, ari))

}


# ==============================================================================
# Scenario 3: Sample Clustering (Manhattan)
# ==============================================================================
# Same row-data setup as Scenario 1 (dat1 reused — already in scope).
# Ward/Median/Centroid excluded: no centroid concept in L1 space.
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Scenario 3: Sample Clustering (Manhattan)\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-22s  ARI\n", "Linkage Method"))
cat(strrep("-", 37), "\n")

for (lm in linkage_methods_gen) {

  linkage_mat <- agglomerative_cluster(
    data     = dat1$X,
    method   = lm$method,
    metric   = DistanceMetric$Manhattan,
    use_mmap = FALSE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = n1, num_clusters = k1)
  ari         <- compute_ari(dat1$true_labels, pred_labels)
  cat(sprintf("%-22s  %.4f\n", lm$name, ari))

}


# ==============================================================================
# Scenario 4: Variable Clustering (Correlation_LSH_Filter)
# ==============================================================================
# Same block-correlated setup as Scenario 2 (X2_t / dat2 reused — in scope).
#
# Exact correlation is computed only when the Hamming distance between 64-bit
# SimHash signatures is <= 14; otherwise the distance is set to 1.0.
# Ward/Median/Centroid excluded (no special projected path for LSH_Filter).
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Scenario 4: Variable Clustering (Correlation_LSH_Filter)\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-22s  ARI\n", "Linkage Method"))
cat(strrep("-", 37), "\n")

for (lm in linkage_methods_gen) {

  linkage_mat <- agglomerative_cluster(
    data     = X2_t,
    method   = lm$method,
    metric   = DistanceMetric$Correlation_LSH_Filter,
    use_mmap = FALSE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = p2, num_clusters = k2)
  ari         <- compute_ari(dat2$true_labels, pred_labels)
  cat(sprintf("%-22s  %.4f\n", lm$name, ari))

}


# ==============================================================================
# Scenario 5: Variable Clustering (Correlation_LSH_Approx)
# ==============================================================================
# Same block-correlated setup as Scenario 2 (X2_t / dat2 reused — in scope).
#
# O(1) distance approximation: 1 - |cos(pi * hamming / 64)|, where hamming is
# the popcount difference between 64-bit SimHash signatures.
#
# Ward is permitted here: the dispatcher routes Ward + Correlation_LSH_Approx
# to ProjectedGeometricUpdatePolicy, which performs Ward clustering in a 64D
# SimHash-projected Euclidean space.
# Note: accuracy may degrade for spatially correlated (e.g. AR(1)) data
# where random-hyperplane quantization shatters tight clusters.
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Scenario 5: Variable Clustering (Correlation_LSH_Approx)\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-22s  ARI\n", "Linkage Method"))
cat(strrep("-", 37), "\n")

for (lm in linkage_methods_lsh_approx) {

  linkage_mat <- agglomerative_cluster(
    data     = X2_t,
    method   = lm$method,
    metric   = DistanceMetric$Correlation_LSH_Approx,
    use_mmap = FALSE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = p2, num_clusters = k2)
  ari         <- compute_ari(dat2$true_labels, pred_labels)
  cat(sprintf("%-22s  %.4f\n", lm$name, ari))

}


cat("\n", strrep("=", 70), "\n", sep = "")
cat("demo_agg_hac_01_in_memory.R complete.\n")
