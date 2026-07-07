# ==============================================================================
# demo_agg_hac_02_mmap.R
# ==============================================================================
#
# Mirrors demo_agg_hac_01_in_memory.R but routes all clustering through
# mmap-backed storage (use_mmap = TRUE).
#
# Scenario 1: Sample Clustering (Euclidean)
#   Same structured row data as demo_01; confirms ARI is identical via mmap path.
#   Euclidean-geometry linkage set: Ward, Average, Complete, Single.
#
# Scenario 2: Variable Clustering (Correlation)
#   Same latent-factor block data as demo_01; t(X) written to mmap.
#   General linkage set: Average, Complete, Single, WPGMA.
#   (Ward/Median/Centroid require Euclidean geometry and are excluded.)
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
cat("Agglomerative Hierarchical Clustering (memory-mapped)\n")
cat(strrep("=", 70), "\n\n")

# Linkage methods: Euclidean-geometry set. Used for Scenario 1 only.
linkage_methods_eucl <- list(
  list(name = "Ward",            method = LinkageMethod$Ward),
  list(name = "Average (UPGMA)", method = LinkageMethod$Average),
  list(name = "Complete",        method = LinkageMethod$Complete),
  list(name = "Single (SLINK)",  method = LinkageMethod$Single)
)

# Linkage methods: general set — purely distance-based, metric-agnostic.
# Used for Scenario 2.
linkage_methods_gen <- list(
  list(name = "Average (UPGMA)", method = LinkageMethod$Average),
  list(name = "Complete",        method = LinkageMethod$Complete),
  list(name = "Single (SLINK)",  method = LinkageMethod$Single),
  list(name = "WPGMA",           method = LinkageMethod$WPGMA)
)


# ==============================================================================
# Scenario 1: Sample Clustering (Euclidean) — mmap path
# ==============================================================================
# Same structured data as demo_01 Scenario 1. X (n1 x p1) is written to a
# temporary mmap file; clustering operates with use_mmap = TRUE.
# Expected: ARI identical to the in-memory result.
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("Scenario 1: Sample Clustering (Euclidean) — mmap\n")
cat(strrep("=", 70), "\n\n")

n1        <- 300L
p1        <- 5000L
k1        <- 3L
centers1  <- c(0.0, 5.0, 20.0)
std_devs1 <- c(1.0, 1.0, 1.0)

mmap_path1 <- tempfile(fileext = ".bin")
on.exit(unlink(mmap_path1), add = TRUE)

cat(sprintf("Generating data: n = %d, p = %d, k = %d ...\n", n1, p1, k1))
dat1   <- generate_sample_clusters(n1, p1, k1, centers1, std_devs1, seed = 42L)
X1_mmap <- convert_to_memory_mapped(dat1$X, mmap_path1)
cat("Done. Data written to mmap file.\n\n")

cat(sprintf("%-20s  ARI\n", "Linkage Method"))
cat(strrep("-", 35), "\n")

for (lm in linkage_methods_eucl) {

  linkage_mat <- agglomerative_cluster(
    data     = X1_mmap,
    method   = lm$method,
    metric   = DistanceMetric$Euclidean,
    use_mmap = TRUE
  )

  pred_labels <- cut_tree(linkage_mat, num_orig_objs = n1, num_clusters = k1)
  ari         <- compute_ari(dat1$true_labels, pred_labels)
  cat(sprintf("%-20s  %.4f\n", lm$name, ari))

}


# ==============================================================================
# Scenario 2: Variable Clustering (Correlation) — mmap path
# ==============================================================================
# Same structured data as demo_01 Scenario 2. The standardized transpose
# t(X) (p2 x n2) is written to a second mmap file so that rows correspond
# to the variables being clustered.
# Expected: ARI identical to the in-memory result.
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Scenario 2: Variable Clustering (Correlation) — mmap\n")
cat(strrep("=", 70), "\n\n")

n2           <- 1000L
p2           <- 5000L
k2           <- 5L
block_sizes2 <- c(2000L, 1000L, 800L, 700L, 500L)

mmap_path2 <- tempfile(fileext = ".bin")
on.exit(unlink(mmap_path2), add = TRUE)

cat(sprintf("Generating data: n = %d, p = %d, k = %d ...\n", n2, p2, k2))
dat2   <- generate_variable_blocks(n2, p2, k2, block_sizes2,
                                   signal_strength = 0.85,
                                   noise_strength  = 0.15,
                                   seed = 42L)
cat("Standardizing columns for correlation ...\n")
X2_t    <- t(standardize_for_correlation(dat2$X))  # p2 x n2
X2_mmap <- convert_to_memory_mapped(X2_t, mmap_path2)
cat("Done. Transposed data written to mmap file.\n\n")

cat(sprintf("%-20s  ARI\n", "Linkage Method"))
cat(strrep("-", 35), "\n")

for (lm in linkage_methods_gen) {

  linkage_mat <- agglomerative_cluster(
    data     = X2_mmap,
    method   = lm$method,
    metric   = DistanceMetric$Correlation,
    use_mmap = TRUE
  )
  pred_labels <- cut_tree(linkage_mat, num_orig_objs = p2, num_clusters = k2)
  ari         <- compute_ari(dat2$true_labels, pred_labels)
  cat(sprintf("%-20s  %.4f\n", lm$name, ari))

}


cat("\n", strrep("=", 70), "\n", sep = "")
cat("demo_agg_hac_02_mmap.R complete.\n")
