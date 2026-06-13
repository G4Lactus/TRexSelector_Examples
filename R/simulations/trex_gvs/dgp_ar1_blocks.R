# ==============================================================================
# dgp_ar1_blocks.R
# ==============================================================================
# Data-generating process for the block-structured AR(1) correlation scenario.
# Sourced by demo_trex_gvs_10.R.
#
# Contents:
#   generate_ar1_blocks(n, p, rho, sigma_y)
#       Interactive demo DGP — 4 blocks (3 active, 1 trap), AR(1) within-block.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_ar1_blocks
# ------------------------------------------------------------------------------
# DGP Scenario: Block-Structured AR(1) Correlation.
# Generates distinct blocks where correlation decays exponentially: Cor(i,j) = rho^|i-j|.
# 4 blocks (sizes 20, 50, 80, 65): 3 active, 1 inactive trap; shuffled with noise gaps.

generate_ar1_blocks <- function(n = 200, p = 500, rho = 0.8, sigma_y = 15) {

  group_sizes <- c(20, 50, 80, 65)
  is_active <- c(TRUE, TRUE, TRUE, FALSE) # Block 4 (size 65) is the inactive trap
  num_blocks <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  if (p < total_grouped + num_blocks - 1) {
    stop("Error: p is too small to hold all blocks and separating noise.")
  }

  # 1. Shuffle block order
  block_order <- sample(1:num_blocks)

  # 2. Distribute separating white noise
  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1)
  gaps <- as.numeric(rmultinom(1, size = rem_noise, prob = rep(1, num_blocks + 1)))
  for(i in 2:num_blocks) gaps[i] <- gaps[i] + 1

  # Initialize
  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)

  current_idx <- 1
  block_indices <- list()

  # 3. Build the matrix sequentially
  for (i in 1:(num_blocks + 1)) {
    # Place noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rnorm(n * gaps[i], mean = 0, sd = 1)
      current_idx <- end_idx + 1
    }

    # Place the AR(1) block
    if (i <= num_blocks) {
      b_id <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1

      # Build the AR(1) Covariance Matrix for this block
      Sigma_ar1 <- matrix(0, nrow = b_size, ncol = b_size)
      for (r in 1:b_size) {
        for (c in 1:b_size) {
          Sigma_ar1[r, c] <- rho^(abs(r - c))
        }
      }

      # Generate multivariate normal data using Cholesky decomposition
      # X_block = Z * chol(Sigma), where Z is n x b_size standard normal
      Z <- matrix(rnorm(n * b_size), nrow = n, ncol = b_size)
      X[, current_idx:end_idx] <- Z %*% chol(Sigma_ar1)

      if (is_active[b_id]) {
        beta_true[current_idx:end_idx] <- 3
      }

      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1
    }
  }

  # 4. Generate response y strictly via linear model
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    block_indices = block_indices
  ))
}


# ==============================================================================
# dgp_ar1_blocks_snr
# ==============================================================================
#' SNR-calibrated block-structured AR(1) DGP.
#'
#' 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
#' Within each block: Cor(X_j, X_k) = rho^|j-k|.
#' Active set: s = 150.
#'
#' @param n   Number of observations.
#' @param p   Total predictors. Default 500.
#' @param snr Signal-to-noise ratio.
#' @param rho AR(1) coefficient (lag-1 correlation). Default 0.85.
#' @param seed RNG seed. Default NULL.
#' @return list(X, y, beta, n, p, s, snr, sigma_y, rho, block_indices, block_order)
dgp_ar1_blocks_snr <- function(n, p = 500, snr, rho = 0.85, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  group_sizes   <- c(20L, 50L, 80L, 65L)
  is_active     <- c(TRUE, TRUE, TRUE, FALSE)
  num_blocks    <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  block_order <- sample(seq_len(num_blocks))

  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1L)
  gaps      <- as.numeric(rmultinom(1, size = rem_noise,
                                    prob = rep(1, num_blocks + 1)))
  for (i in 2:num_blocks) gaps[i] <- gaps[i] + 1L

  X             <- matrix(0, nrow = n, ncol = p)
  beta          <- rep(0, p)
  current_idx   <- 1L
  block_indices <- vector("list", num_blocks)

  for (i in seq_len(num_blocks + 1L)) {
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1L
      X[, current_idx:end_idx] <- rnorm(n * gaps[i])
      current_idx <- end_idx + 1L
    }
    if (i <= num_blocks) {
      b_id   <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1L

      # AR(1) covariance: Sigma[r,c] = rho^|r-c|
      r_idx    <- seq_len(b_size)
      Sigma_ar1 <- rho^abs(outer(r_idx, r_idx, `-`))
      Z         <- matrix(rnorm(n * b_size), nrow = n, ncol = b_size)
      X[, current_idx:end_idx] <- Z %*% chol(Sigma_ar1)

      if (is_active[b_id]) beta[current_idx:end_idx] <- 3
      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1L
    }
  }

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, sd = sigma_y)

  list(
    X             = X,
    y             = y,
    beta          = beta,
    n             = n,
    p             = p,
    s             = 150L,
    snr           = snr,
    sigma_y       = sigma_y,
    rho           = rho,
    block_indices = block_indices,
    block_order   = block_order
  )
}
