# ==============================================================================
# hac_helpers.R
# ==============================================================================
#
# Shared helper functions for agglomerative hierarchical clustering demos.
# Sourced by demo_agg_hac_01_in_memory.R and demo_agg_hac_02_mmap.R.
#
# Functions:
#   - compute_ari()              Adjusted Rand Index (no external dependencies)
#   - generate_sample_clusters() Structured row data with distinct Gaussian centers
#   - generate_variable_blocks() Latent-factor block-correlated column data
#   - standardize_for_correlation() Mean-center + L2-normalize each column
#
# ==============================================================================


# ==============================================================================
# compute_ari
# ==============================================================================

#' @title Adjusted Rand Index
#'
#' @description Computes the Adjusted Rand Index (ARI) between two label vectors.
#'   R port of the inline C++ implementation in demo_mlm_hac_01.cpp.
#'
#' @param true_labels Integer vector of ground-truth cluster labels.
#' @param pred_labels Integer vector of predicted cluster labels (same length).
#'
#' @return ARI score in [-1, 1]; 1.0 = perfect agreement.
compute_ari <- function(true_labels, pred_labels) {

  n <- length(true_labels)
  if (n == 0L || n != length(pred_labels)) return(0.0)

  choose2 <- function(x) x * (x - 1.0) / 2.0

  tbl        <- table(true = true_labels, pred = pred_labels)
  sum_comb_c <- sum(choose2(as.double(tbl)))
  sum_comb_a <- sum(choose2(as.double(rowSums(tbl))))
  sum_comb_b <- sum(choose2(as.double(colSums(tbl))))

  n_choose_2 <- choose2(as.double(n))
  expected   <- (sum_comb_a * sum_comb_b) / n_choose_2
  max_index  <- (sum_comb_a + sum_comb_b) / 2.0

  if (max_index == expected) return(1.0)
  (sum_comb_c - expected) / (max_index - expected)
}


# ==============================================================================
# generate_sample_clusters
# ==============================================================================

#' @title Generate row-clustered data with distinct Gaussian centers
#'
#' @description Generates an (n x p) matrix where rows are partitioned into k
#'   equal-size clusters. Each cluster has a distinct location-scale transform
#'   applied: X_{i,j} = center_k + std_dev_k * Z, Z ~ N(0,1).
#'   R port of generate_hierarchical_row_data() in demo_mlm_hac_01.cpp.
#'
#' @param n       Total number of samples. Must be divisible by k.
#' @param p       Number of features per sample.
#' @param k       Number of clusters.
#' @param centers Numeric vector of length k. Mean offset per cluster.
#' @param std_devs Numeric vector of length k. Standard deviation per cluster.
#' @param seed    Random seed for reproducibility.
#'
#' @return A list with:
#'   \item{X}{Numeric matrix (n x p).}
#'   \item{true_labels}{Integer vector of length n; 0-based cluster labels.}
generate_sample_clusters <- function(n, p, k, centers, std_devs, seed = 42L) {

  if (n %% k != 0L)           stop("n must be divisible by k")
  if (length(centers)  != k)  stop("centers must have length k")
  if (length(std_devs) != k)  stop("std_devs must have length k")

  cluster_size <- n %/% k
  true_labels  <- rep(seq_len(k) - 1L, each = cluster_size)  # 0-based

  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  for (ki in seq_len(k)) {
    rows     <- ((ki - 1L) * cluster_size + 1L):(ki * cluster_size)
    X[rows, ] <- centers[ki] + std_devs[ki] * X[rows, ]
  }

  list(X = X, true_labels = true_labels)
}


# ==============================================================================
# generate_variable_blocks
# ==============================================================================

#' @title Generate block-correlated variables via a latent factor model
#'
#' @description Generates an (n x p) matrix where columns (variables) are
#'   partitioned into k blocks. Within each block, all columns share a latent
#'   signal drawn from N(0,1), producing block-diagonal Pearson correlation.
#'   X_{i,j} = signal_strength * latent_{i,k} + noise_strength * eps_{i,j}
#'   R port of generate_correlated_variables() in demo_mlm_hac_01.cpp.
#'
#' @param n              Number of observations (rows).
#' @param p              Total number of variables (columns). Must equal sum(block_sizes).
#' @param k              Number of variable blocks.
#' @param block_sizes    Integer vector of length k. Number of variables per block.
#' @param signal_strength Scalar in (0,1). Shared-signal weight. Default 0.85.
#' @param noise_strength  Scalar in (0,1). Independent-noise weight. Default 0.15.
#' @param seed           Random seed for reproducibility.
#'
#' @return A list with:
#'   \item{X}{Numeric matrix (n x p).}
#'   \item{true_labels}{Integer vector of length p; 0-based block labels for columns.}
generate_variable_blocks <- function(n, p, k, block_sizes,
                                     signal_strength = 0.85,
                                     noise_strength  = 0.15,
                                     seed = 42L) {
  if (length(block_sizes) != k)  stop("block_sizes must have length k")
  if (sum(block_sizes)    != p)  stop("sum(block_sizes) must equal p")

  true_labels <- rep(seq_len(k) - 1L, times = block_sizes)  # 0-based, length p

  set.seed(seed)
  latent_signals <- matrix(rnorm(n * k), nrow = n, ncol = k)  # n x k

  X         <- matrix(0.0, nrow = n, ncol = p)
  col_start <- 1L
  for (ki in seq_len(k)) {
    col_end   <- col_start + block_sizes[ki] - 1L
    noise_blk <- matrix(rnorm(n * block_sizes[ki]), nrow = n, ncol = block_sizes[ki])
    X[, col_start:col_end] <- signal_strength * latent_signals[, ki] + noise_strength  * noise_blk
    col_start <- col_end + 1L
  }

  list(X = X, true_labels = true_labels)
}


# ==============================================================================
# standardize_for_correlation
# ==============================================================================

#' @title Mean-center and L2-normalize each column
#'
#' @description Transforms each column of X so that it has zero mean and unit
#'   L2 norm. Required before Pearson correlation clustering.
#'   R port of standardize_for_correlation() in demo_mlm_hac_01.cpp.
#'
#' @param X Numeric matrix (n x p).
#'
#' @return The standardized matrix (n x p).
standardize_for_correlation <- function(X) {
  X         <- sweep(X, 2L, colMeans(X), "-")
  col_norms <- sqrt(colSums(X^2))
  col_norms[col_norms < 1e-12] <- 1.0  # guard against zero-variance columns
  sweep(X, 2L, col_norms, "/")
}
