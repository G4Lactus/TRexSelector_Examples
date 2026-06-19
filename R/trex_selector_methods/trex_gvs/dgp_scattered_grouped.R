# ==============================================================================
# dgp_scattered_grouped.R
# ==============================================================================
# Data-generating process for the scattered-grouped predictors scenario.
# Sourced by demo_trex_gvs_03_hastie_realistic.R.
#
# Contents:
#   generate_scattered_grouped_data(n, p, n_groups, group_size, sigma_y)
#       Interactive demo DGP — groups randomly distributed across column indices.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_scattered_grouped_data
# ------------------------------------------------------------------------------
# Realistic DGP: Scattered Grouped Variables for Elastic Net.
# Groups are randomly distributed across the column indices.

generate_scattered_grouped_data <- function(n = 200, p = 500, n_groups = 3, group_size = 50, sigma_y = 15) {

  if (n_groups * group_size > p) {
    stop("Error: Total number of grouped variables exceeds p.")
  }

  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)

  # 1. Randomly partition the indices into groups and noise
  all_indices <- 1:p
  group_indices <- list()

  for (k in 1:n_groups) {
    # Sample indices for group k and remove them from the available pool
    idx <- sample(all_indices, size = group_size, replace = FALSE)
    group_indices[[k]] <- sort(idx)
    all_indices <- setdiff(all_indices, idx)
  }

  # The remaining indices are assigned as pure noise
  noise_indices <- sort(all_indices)

  # 2. Generate the data
  sd_x <- sqrt(0.01) # within-group noise (correlation ~ 0.99)

  # Generate scattered grouped predictors
  for (k in 1:n_groups) {
    Z_k <- rnorm(n, mean = 0, sd = 1) # Latent factor for group k

    for (i in group_indices[[k]]) {
      X[, i] <- Z_k + rnorm(n, mean = 0, sd = sd_x)
      beta_true[i] <- 3  # All members of the correlated group are active
    }
  }

  # Generate pure noise features
  if (length(noise_indices) > 0) {
    for (i in noise_indices) {
      X[, i] <- rnorm(n, mean = 0, sd = 1)
      beta_true[i] <- 0  # Noise features are null
    }
  }

  # 3. Generate the response variable y
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  # Return data and the index mapping for evaluation
  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    group_indices = group_indices,
    noise_indices = noise_indices
  ))
}


# ------------------------------------------------------------------------------
# dgp_scattered_grouped_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP for the scattered-grouped scenario with SNR-calibrated noise.
#
# n_groups groups of group_size variables each, randomly distributed across p
# columns, sharing a latent factor within each group (within-group correlation
# rho = 1 / (1 + sd_x^2)).
#
# The noise standard deviation sigma_y is derived from the target SNR:
#
#   signal   = X %*% beta
#   sigma_y  = sqrt( var(signal) / snr )      (n-1 denominator)
#
# @param n          Number of observations.
# @param p          Number of predictors. Default 500.
# @param snr        Target signal-to-noise ratio Var(signal) / sigma_y^2.
# @param sd_x       Within-group noise sd. Controls rho = 1/(1+sd_x^2).
#                   Default sqrt(0.01) => rho ~ 0.99.
# @param n_groups   Number of active groups. Default 3.
# @param group_size Number of variables per group. Default 50.
# @param seed       Integer RNG seed. Default NULL.
#
# @return Named list:
#   X           n x p predictor matrix
#   y           length-n response
#   beta        length-p true coefficient vector (field consistent with .run_mc)
#   n, p, s     dimensions (s = n_groups * group_size)
#   snr         target SNR
#   sigma_y     realised noise standard deviation
#   sd_x        within-group noise sd used
#   n_groups    number of active groups
#   group_size  variables per group

dgp_scattered_grouped_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                                      n_groups = 3L, group_size = 50L,
                                      seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (n_groups * group_size > p) {
    stop("'n_groups * group_size' must not exceed 'p'.")
  }

  # Randomly partition indices into groups and noise
  all_indices   <- seq_len(p)
  group_indices <- vector("list", n_groups)

  for (k in seq_len(n_groups)) {
    idx              <- sample(all_indices, size = group_size, replace = FALSE)
    group_indices[[k]] <- sort(idx)
    all_indices      <- setdiff(all_indices, idx)
  }

  # True coefficient vector — all group members are active (beta = 3)
  beta <- rep(0, p)
  X    <- matrix(0, nrow = n, ncol = p)

  for (k in seq_len(n_groups)) {
    Z_k <- rnorm(n, mean = 0, sd = 1)
    for (i in group_indices[[k]]) {
      X[, i]   <- Z_k + rnorm(n, mean = 0, sd = sd_x)
      beta[i]  <- 3
    }
  }

  # Remaining columns: independent white noise
  noise_indices <- sort(all_indices)
  if (length(noise_indices) > 0) {
    for (i in noise_indices) X[, i] <- rnorm(n, mean = 0, sd = 1)
  }

  # SNR-calibrated noise
  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)

  y <- signal + rnorm(n, mean = 0, sd = sigma_y)

  list(
    X          = X,
    y          = y,
    beta       = beta,
    n          = n,
    p          = p,
    s          = as.integer(n_groups * group_size),
    snr        = snr,
    sigma_y    = sigma_y,
    sd_x       = sd_x,
    n_groups   = n_groups,
    group_size = group_size
  )
}
