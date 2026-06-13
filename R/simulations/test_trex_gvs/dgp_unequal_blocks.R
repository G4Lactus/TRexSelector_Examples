# ==============================================================================
# dgp_unequal_blocks.R
# ==============================================================================
# Data-generating process for the unequal-blockwise predictors scenario.
# Sourced by demo_trex_gvs_04_hastie_realistic_2.R.
#
# Contents:
#   generate_unequal_blockwise_data(n, p, sigma_y)
#       Interactive demo DGP — three contiguous unequal blocks (sizes 20, 50, 80).
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_unequal_blockwise_data
# ------------------------------------------------------------------------------
# DGP Scenario: Unequal Blockwise Groups.
# Three groups of sizes 20, 50, and 80 (Sum = 150) placed contiguously.
# All 150 variables in these groups are active (beta = 3).
# The remainder are pure white noise (beta = 0).

generate_unequal_blockwise_data <- function(n = 200, p = 500, sigma_y = 15) {

  # Define the unequal group sizes
  group_sizes <- c(20, 50, 80)
  total_grouped <- sum(group_sizes)

  if (p < total_grouped) {
    stop("Error: p must be at least large enough to hold all groups.")
  }

  # 1. Define the true beta coefficient vector
  # First 150 are active, the rest are 0
  beta_true <- rep(0, p)
  beta_true[1:total_grouped] <- 3

  X <- matrix(0, nrow = n, ncol = p)
  sd_x <- sqrt(0.01) # within-group noise

  # 2. Generate the unequal blocks
  current_idx <- 1

  for (k in 1:length(group_sizes)) {
    size_k <- group_sizes[k]
    Z_k <- rnorm(n, mean = 0, sd = 1) # Latent factor for group k

    # Fill the block
    end_idx <- current_idx + size_k - 1
    for (i in current_idx:end_idx) {
      X[, i] <- Z_k + rnorm(n, mean = 0, sd = sd_x)
    }

    current_idx <- end_idx + 1
  }

  # 3. Generate pure noise for the remainder
  if (p > total_grouped) {
    for (i in current_idx:p) {
      X[, i] <- rnorm(n, mean = 0, sd = 1)
    }
  }

  # 4. Generate the response variable y
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  # Return data
  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    group_sizes = group_sizes
  ))
}


# ------------------------------------------------------------------------------
# dgp_unequal_blocks_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP: three contiguous unequal blocks (sizes 20, 50, 80).
# All 150 variables are active (beta = 3). Within-block correlation
# rho = 1 / (1 + sd_x^2).  sigma_y is derived from the target SNR:
#
#   signal  = X %*% beta
#   sigma_y = sqrt( var(signal) / snr )    (n-1 denominator)
#
# @param n      Number of observations.
# @param p      Number of predictors. Default 500.
# @param snr    Target signal-to-noise ratio.
# @param sd_x   Within-block noise sd. Default sqrt(0.01) => rho ~ 0.99.
# @param seed   Integer RNG seed. Default NULL.
#
# @return Named list: X, y, beta, n, p, s=150, snr, sigma_y, sd_x, group_sizes.

dgp_unequal_blocks_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                                   seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  group_sizes <- c(20L, 50L, 80L)
  s           <- sum(group_sizes)

  beta         <- rep(0, p)
  beta[1:s]    <- 3

  X           <- matrix(0, nrow = n, ncol = p)
  current_idx <- 1L
  for (k in seq_along(group_sizes)) {
    Z_k     <- rnorm(n, 0, 1)
    end_idx <- current_idx + group_sizes[k] - 1L
    for (i in current_idx:end_idx) X[, i] <- Z_k + rnorm(n, 0, sd_x)
    current_idx <- end_idx + 1L
  }
  for (i in seq(current_idx, p)) X[, i] <- rnorm(n, 0, 1)

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, 0, sigma_y)

  list(
    X           = X,
    y           = y,
    beta        = beta,
    n           = n,
    p           = p,
    s           = as.integer(s),
    snr         = snr,
    sigma_y     = sigma_y,
    sd_x        = sd_x,
    group_sizes = group_sizes
  )
}
