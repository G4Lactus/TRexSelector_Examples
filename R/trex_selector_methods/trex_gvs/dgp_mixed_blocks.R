# ==============================================================================
# dgp_mixed_blocks.R
# ==============================================================================
# Data-generating process for the mixed active/inactive blocks scenario.
# Sourced by demo_trex_gvs_03.R and simulation_trex_gvs_03.R.
#
# Contents:
#   generate_mixed_blocks_data(n, p, sigma_y)
#       Interactive demo DGP — 3 active + 1 inactive block, randomly ordered.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_mixed_blocks_data
# ------------------------------------------------------------------------------
# DGP Scenario: Active and Inactive Correlated Blocks.
# 4 Blocks total.
# Active Blocks: Sizes 20, 50, 80 (beta = 3)
# Inactive Block: Size 65 (beta = 0)
# Randomly ordered, separated by white noise.

generate_mixed_blocks_data <- function(n = 200, p = 500, sigma_y = 15) {

  # Block definitions
  block_sizes <- c(20, 50, 80, 65)
  is_active <- c(TRUE, TRUE, TRUE, FALSE) # Block 4 (size 65) is inactive
  num_blocks <- length(block_sizes)
  total_grouped <- sum(block_sizes)

  if (p < total_grouped + num_blocks - 1) {
    stop("Error: p is too small to hold all blocks and separating noise.")
  }

  # 1. Shuffle the order of the blocks
  block_order <- sample(1:num_blocks)

  # 2. Distribute the noise into (num_blocks + 1) gaps
  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_blocks - 1)

  # Allocate remaining noise across 5 gaps, adding 1 to internal gaps to ensure separation
  gaps <- as.numeric(rmultinom(1, size = rem_noise, prob = rep(1, num_blocks + 1)))
  for(i in 2:num_blocks) {
    gaps[i] <- gaps[i] + 1
  }

  # Initialize
  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)
  sd_x <- sqrt(0.01)

  # Latent driving factors for all 4 blocks
  latent_factors <- list(
    rnorm(n, mean = 0, sd = 1), # Block 1 (Active, 20)
    rnorm(n, mean = 0, sd = 1), # Block 2 (Active, 50)
    rnorm(n, mean = 0, sd = 1), # Block 3 (Active, 80)
    rnorm(n, mean = 0, sd = 1)  # Block 4 (Inactive, 65)
  )

  # 3. Build the matrix sequentially
  current_idx <- 1
  block_indices <- list()

  for (i in 1:(num_blocks + 1)) {
    # Place noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rnorm(n * gaps[i], mean = 0, sd = 1)
      current_idx <- end_idx + 1
    }

    # Place the block
    if (i <= num_blocks) {
      b_id <- block_order[i]
      b_size <- block_sizes[b_id]
      end_idx <- current_idx + b_size - 1

      Z <- latent_factors[[b_id]]
      for (c in current_idx:end_idx) {
        X[, c] <- Z + rnorm(n, mean = 0, sd = sd_x)
      }

      # Apply beta_true ONLY if the block is designated as active
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


# ------------------------------------------------------------------------------
# dgp_mixed_blocks_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP: 3 active blocks (20, 50, 80) + 1 inactive block (65),
# randomly ordered and separated by white noise.
# Active set: s = 150 (blocks 1-3 only). Within-block correlation
# rho = 1 / (1 + sd_x^2).  sigma_y derived from target SNR.
#
# @param n     Number of observations.
# @param p     Number of predictors. Default 500.
# @param snr   Target signal-to-noise ratio.
# @param sd_x  Within-block noise sd. Default sqrt(0.01) => rho ~ 0.99.
# @param seed  Integer RNG seed. Default NULL.
#
# @return Named list: X, y, beta, n, p, s=150, snr, sigma_y, sd_x,
#                     block_sizes, is_active, block_order, block_indices.

dgp_mixed_blocks_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  block_sizes   <- c(20L, 50L, 80L, 65L)
  is_active     <- c(TRUE, TRUE, TRUE, FALSE)
  num_blocks    <- length(block_sizes)
  total_grouped <- sum(block_sizes)

  block_order <- sample(seq_len(num_blocks))
  num_noise   <- p - total_grouped
  rem_noise   <- num_noise - (num_blocks - 1L)
  gaps        <- as.numeric(rmultinom(1, size = rem_noise,
                                      prob = rep(1, num_blocks + 1)))
  for (i in 2:num_blocks) gaps[i] <- gaps[i] + 1

  X             <- matrix(0, nrow = n, ncol = p)
  beta          <- rep(0, p)
  latent        <- lapply(seq_len(num_blocks), function(k) rnorm(n, 0, 1))
  current_idx   <- 1L
  block_indices <- vector("list", num_blocks)

  for (i in 1:(num_blocks + 1L)) {
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1L
      X[, current_idx:end_idx] <- matrix(rnorm(n * gaps[i], 0, 1), n, gaps[i])
      current_idx <- end_idx + 1L
    }
    if (i <= num_blocks) {
      b_id    <- block_order[i]
      b_size  <- block_sizes[b_id]
      end_idx <- current_idx + b_size - 1L
      Z       <- latent[[b_id]]
      for (c in current_idx:end_idx) X[, c] <- Z + rnorm(n, 0, sd_x)
      if (is_active[b_id]) beta[current_idx:end_idx] <- 3
      block_indices[[b_id]] <- current_idx:end_idx
      current_idx <- end_idx + 1L
    }
  }

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, 0, sigma_y)

  list(
    X             = X,
    y             = y,
    beta          = beta,
    n             = n,
    p             = p,
    s             = as.integer(sum(block_sizes[is_active])),
    snr           = snr,
    sigma_y       = sigma_y,
    sd_x          = sd_x,
    block_sizes   = block_sizes,
    is_active     = is_active,
    block_order   = block_order,
    block_indices = block_indices
  )
}
