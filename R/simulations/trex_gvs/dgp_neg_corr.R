# ==============================================================================
# dgp_neg_corr.R
# ==============================================================================
# Data-generating process for the strict linear negative correlations scenario.
# Sourced by demo_trex_gvs_06_hastie_2.R.
#
# Contents:
#   generate_strict_linear_negative_correlations(n, p, sigma_y)
#       Interactive demo DGP — active group with +Z/-Z structure; inactive trap.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_strict_linear_negative_correlations
# ------------------------------------------------------------------------------
# DGP Scenario: Strictly Scaled Idealized Example (Linear Model Version).
# Predictors contain negative correlations (+Z and -Z).
# Group 1 (indices 1-100): Active (+Z1/-Z1, beta=+3/-3).
# Group 2 (indices 101-200): Inactive Trap (+Z2/-Z2, beta=0).
# Indices 201-p: white noise.

generate_strict_linear_negative_correlations <- function(n = 200, p = 500, sigma_y = 15) {

  # 1. Generate independent latent variables
  Z1 <- rnorm(n, mean = 0, sd = 1)
  Z2 <- rnorm(n, mean = 0, sd = 1)

  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)
  sd_x <- sqrt(0.01)

  g_size <- 100

  # 2. Build Group 1 (The Active Group) - Indices 1 to 100
  # First half: +Z1
  for (i in 1:(g_size/2)) {
    X[, i] <- Z1 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 3
  }
  # Second half: -Z1 (Negatively correlated with first half)
  for (i in (g_size/2 + 1):g_size) {
    X[, i] <- -Z1 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- -3 # Negative beta perfectly counteracts the negative X
  }

  # 3. Build Group 2 (The Inactive Trap Group) - Indices 101 to 200
  # First half: +Z2
  for (i in 1:(g_size/2)) {
    X[, g_size + i] <- Z2 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[g_size + i] <- 0
  }
  # Second half: -Z2
  for (i in (g_size/2 + 1):g_size) {
    X[, g_size + i] <- -Z2 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[g_size + i] <- 0
  }

  # 4. Build High-Dimensional Noise - Indices 201 to p
  for (i in (2*g_size + 1):p) {
    X[, i] <- rnorm(n, mean = 0, sd = 1)
    beta_true[i] <- 0
  }

  # 5. Generate response y STRICTLY via linear model
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X = X,
    y = y,
    beta_true = beta_true
  ))
}


# ------------------------------------------------------------------------------
# dgp_neg_corr_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP: strictly scaled negative-correlation scenario.
# Group 1 (1-100): active, +Z1/-Z1, beta = +3/-3.  s = 100.
# Group 2 (101-200): inactive trap, +Z2/-Z2, beta = 0.
# Indices 201-p: white noise.
# Within-half-group correlation rho = 1/(1+sd_x^2); cross-half rho ~ -rho.
# sigma_y derived from target SNR.
#
# @param n     Number of observations.
# @param p     Number of predictors. Default 500.
# @param snr   Target signal-to-noise ratio.
# @param sd_x  Within-group noise sd. Default sqrt(0.01) => rho ~ 0.99.
# @param seed  Integer RNG seed. Default NULL.
#
# @return Named list: X, y, beta, n, p, s=100, snr, sigma_y, sd_x.

dgp_neg_corr_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  Z1 <- rnorm(n, 0, 1)
  Z2 <- rnorm(n, 0, 1)

  X      <- matrix(0, nrow = n, ncol = p)
  beta   <- rep(0, p)
  g_size <- 100L

  # Group 1: active (+Z1/-Z1)
  for (i in seq_len(g_size / 2L)) {
    X[, i] <- Z1 + rnorm(n, 0, sd_x); beta[i] <- 3
  }
  for (i in seq(g_size / 2L + 1L, g_size)) {
    X[, i] <- -Z1 + rnorm(n, 0, sd_x); beta[i] <- -3
  }

  # Group 2: inactive trap (+Z2/-Z2)
  for (i in seq_len(g_size / 2L)) {
    X[, g_size + i] <- Z2 + rnorm(n, 0, sd_x)
  }
  for (i in seq(g_size / 2L + 1L, g_size)) {
    X[, g_size + i] <- -Z2 + rnorm(n, 0, sd_x)
  }

  # White noise
  if (p > 2L * g_size) {
    for (i in seq(2L * g_size + 1L, p)) X[, i] <- rnorm(n, 0, 1)
  }

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, 0, sigma_y)

  list(
    X       = X,
    y       = y,
    beta    = beta,
    n       = n,
    p       = p,
    s       = as.integer(g_size),
    snr     = snr,
    sigma_y = sigma_y,
    sd_x    = sd_x
  )
}
