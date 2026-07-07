# ==============================================================================
# dgp_t3_equi_blocks.R
# ==============================================================================
# Data-generating process for the heavy-tailed equi-correlated blocks scenario.
# Sourced by demo_trex_gvs_06.R and simulation_trex_gvs_06.R.
#
# Contents:
#   generate_t3_equicorrelated_blocks(n, p, rho, sigma_y, df)
#       Interactive demo DGP — 4 blocks with Student-t distributed predictors.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_t3_equicorrelated_blocks
# ------------------------------------------------------------------------------
# DGP Scenario: Heavy-Tailed Equi-Correlated Blocks.
# Uses Student's t-distribution with df degrees of freedom instead of Gaussian.
# Variables are scaled to var=1 to preserve the correlation 'rho' and SNR.
# 4 blocks (sizes 20, 50, 80, 65): 3 active, 1 inactive trap; shuffled with noise gaps.

generate_t3_equicorrelated_blocks <- function(n = 200, p = 500, rho = 0.8, sigma_y = 15, df = 3) {

  # Helper function to generate scaled t-distributed variables (Variance = 1)
  rt_scaled <- function(n_samples) {
    return(rt(n_samples, df = df) / sqrt(df / (df - 2)))
  }

  group_sizes <- c(20, 50, 80, 65)
  is_active <- c(TRUE, TRUE, TRUE, FALSE) # Block 4 (size 65) is the trap
  num_blocks <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  if (p < total_grouped + num_blocks - 1) {
    stop("Error: p is too small to hold all blocks and separating noise.")
  }

  # 1. Shuffle block order
  block_order <- sample(1:num_blocks)

  # 2. Distribute separating noise (heavy-tailed white noise)
  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1)

  gaps <- as.numeric(rmultinom(1, size = rem_noise, prob = rep(1, num_blocks + 1)))
  for(i in 2:num_blocks) gaps[i] <- gaps[i] + 1

  # Initialize
  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)

  # Latent driving factors for all blocks (Heavy-tailed)
  latent_factors <- list()
  for (i in 1:num_blocks) latent_factors[[i]] <- rt_scaled(n)

  # 3. Build the matrix sequentially
  current_idx <- 1
  block_indices <- list()

  w_Z <- sqrt(rho)
  w_E <- sqrt(1 - rho)

  for (i in 1:(num_blocks + 1)) {
    # Place heavy-tailed noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rt_scaled(n * gaps[i])
      current_idx <- end_idx + 1
    }

    # Place the heavy-tailed equi-correlated block
    if (i <= num_blocks) {
      b_id <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1

      Z <- latent_factors[[b_id]]
      for (c in current_idx:end_idx) {
        # Construct equi-correlated variable with heavy tails
        X[, c] <- w_Z * Z + w_E * rt_scaled(n)
      }

      if (is_active[b_id]) {
        beta_true[current_idx:end_idx] <- 3
      }

      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1
    }
  }

  # 4. Generate response y with heavy-tailed shocks
  epsilon <- sigma_y * rt_scaled(n)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    block_indices = block_indices
  ))
}


# ==============================================================================
# dgp_t3_equi_blocks_snr
# ==============================================================================
#' SNR-calibrated heavy-tailed equi-correlated blocks DGP.
#'
#' Identical block structure to dgp_equi_blocks_snr but with Student-t(df)
#' distributed latent factors and noise (scaled to variance 1).
#' 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
#' Active set: s = 150.
#'
#' @param n   Number of observations.
#' @param p   Total predictors. Default 500.
#' @param snr Signal-to-noise ratio.
#' @param rho Within-block equi-correlation. Default 0.75.
#' @param df  Degrees of freedom for the t distribution. Default 3.
#' @param seed RNG seed. Default NULL.
#' @return list(X, y, beta, n, p, s, snr, sigma_y, rho, df, block_indices, block_order)
dgp_t3_equi_blocks_snr <- function(n, p = 500, snr, rho = 0.75, df = 3,
                                    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  rt_scaled <- function(n_samp) rt(n_samp, df = df) / sqrt(df / (df - 2))

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
  latent_factors <- lapply(seq_len(num_blocks), function(i) rt_scaled(n))

  current_idx   <- 1L
  block_indices <- vector("list", num_blocks)
  w_Z           <- sqrt(rho)
  w_E           <- sqrt(1 - rho)

  for (i in seq_len(num_blocks + 1L)) {
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1L
      X[, current_idx:end_idx] <- rt_scaled(n * gaps[i])
      current_idx <- end_idx + 1L
    }
    if (i <= num_blocks) {
      b_id   <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1L
      Z <- latent_factors[[b_id]]
      for (col in current_idx:end_idx) {
        X[, col] <- w_Z * Z + w_E * rt_scaled(n)
      }
      if (is_active[b_id]) beta[current_idx:end_idx] <- 3
      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1L
    }
  }

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rt_scaled(n) * sigma_y

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
    df            = df,
    block_indices = block_indices,
    block_order   = block_order
  )
}
