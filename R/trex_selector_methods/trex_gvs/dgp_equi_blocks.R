# ==============================================================================
# dgp_equi_blocks.R
# ==============================================================================
# Data-generating process for the equi-correlated block variables scenario.
# Sourced by demo_trex_gvs_05.R and simulation_trex_gvs_05.R.
#
# Contents:
#   generate_equicorrelated_blocks(n, p, rho, sigma_y)
#       Interactive demo DGP — 4 blocks (3 active, 1 trap), equi-correlated, shuffled.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_equicorrelated_blocks
# ------------------------------------------------------------------------------
# DGP Scenario: Equi-Correlated Block Variables.
# Allows you to specify the exact pairwise correlation 'rho' within groups.
# 4 blocks (sizes 20, 50, 80, 65): 3 active, 1 inactive trap; shuffled with noise gaps.

generate_equicorrelated_blocks <- function(n = 200, p = 500, rho = 0.8, sigma_y = 15) {

  group_sizes <- c(20, 50, 80, 65)
  is_active <- c(TRUE, TRUE, TRUE, FALSE) # Block 4 (size 65) is the inactive trap
  num_blocks <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  if (p < total_grouped + num_blocks - 1) {
    stop("Error: p is too small to hold all blocks and separating noise.")
  }

  # 1. Shuffle the order of the blocks
  block_order <- sample(1:num_blocks)

  # 2. Distribute the separating white noise
  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1)

  gaps <- as.numeric(rmultinom(1, size = rem_noise, prob = rep(1, num_blocks + 1)))
  for(i in 2:num_blocks) gaps[i] <- gaps[i] + 1 # Ensure at least 1 noise var between blocks

  # Initialize
  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)

  # Latent driving factors for all blocks
  latent_factors <- list()
  for (i in 1:num_blocks) latent_factors[[i]] <- rnorm(n, mean = 0, sd = 1)

  # 3. Build the matrix sequentially
  current_idx <- 1
  block_indices <- list()

  # Equi-correlation weights
  w_Z <- sqrt(rho)
  w_E <- sqrt(1 - rho)

  for (i in 1:(num_blocks + 1)) {
    # Place noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rnorm(n * gaps[i], mean = 0, sd = 1)
      current_idx <- end_idx + 1
    }

    # Place the equi-correlated block
    if (i <= num_blocks) {
      b_id <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1

      Z <- latent_factors[[b_id]]
      for (c in current_idx:end_idx) {
        # Construct equi-correlated variable
        X[, c] <- w_Z * Z + w_E * rnorm(n, mean = 0, sd = 1)
      }

      if (is_active[b_id]) {
        beta_true[current_idx:end_idx] <- 3
      }

      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1
    }
  }

  # 4. Generate response y
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    block_indices = block_indices,
    block_order = block_order
  ))
}


# ==============================================================================
# dgp_equi_blocks_snr
# ==============================================================================
#' SNR-calibrated equi-correlated blocks DGP.
#'
#' 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
#' Shuffled into p columns with white-noise gaps.
#' Within-block equi-correlation: Cor(X_j, X_k) = rho exactly.
#' Active set: s = 150.
#'
#' @param n       Number of observations.
#' @param p       Total predictors. Default 500.
#' @param snr     Signal-to-noise ratio (var(signal) / var(noise)).
#' @param rho     Within-block equi-correlation. Default 0.75.
#' @param seed    RNG seed. Default NULL.
#' @return list(X, y, beta, n, p, s, snr, sigma_y, rho, block_indices, block_order)
dgp_equi_blocks_snr <- function(n, p = 500, snr, rho = 0.75, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  group_sizes  <- c(20L, 50L, 80L, 65L)
  is_active    <- c(TRUE, TRUE, TRUE, FALSE)
  num_blocks   <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  block_order <- sample(seq_len(num_blocks))

  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1L)
  gaps      <- as.numeric(rmultinom(1, size = rem_noise,
                                    prob = rep(1, num_blocks + 1)))
  for (i in 2:num_blocks) gaps[i] <- gaps[i] + 1L

  X             <- matrix(0, nrow = n, ncol = p)
  beta          <- rep(0, p)
  latent_factors <- lapply(seq_len(num_blocks), function(i) rnorm(n))

  current_idx   <- 1L
  block_indices <- vector("list", num_blocks)
  w_Z           <- sqrt(rho)
  w_E           <- sqrt(1 - rho)

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
      Z <- latent_factors[[b_id]]
      for (col in current_idx:end_idx) {
        X[, col] <- w_Z * Z + w_E * rnorm(n)
      }
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
