# ==============================================================================
# dgp_neg_traps.R
# ==============================================================================
# Data-generating process for the hardened negative traps scenario.
# Sourced by demo_trex_gvs_07_hastie_3.R.
#
# Contents:
#   generate_hardened_negative_traps(n, p, sigma_y)
#       Interactive demo DGP — active group + 2 inactive +/-Z traps.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_hardened_negative_traps
# ------------------------------------------------------------------------------
# DGP Scenario: Multiple Inactive Traps with Negative Correlations.
# Group 1 (1-100):   Active (+Z1/-Z1, beta=+3/-3).
# Group 2 (101-200): Inactive Trap 1 (+Z2/-Z2, beta=0).
# Noise 1 (201-300): White noise.
# Group 3 (301-360): Inactive Trap 2 (+Z3/-Z3, beta=0).
# Noise 2 (361-p):   White noise.

generate_hardened_negative_traps <- function(n = 200, p = 500, sigma_y = 15) {

  Z1 <- rnorm(n, mean = 0, sd = 1)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  Z3 <- rnorm(n, mean = 0, sd = 1)

  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)
  sd_x <- sqrt(0.01)

  # 1. Group 1: The Active Group (1 to 100)
  for (i in 1:50) {
    X[, i] <- Z1 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 3
  }
  for (i in 51:100) {
    X[, i] <- -Z1 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- -3
  }

  # 2. Group 2: Inactive Trap 1 (101 to 200)
  for (i in 101:150) {
    X[, i] <- Z2 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 0
  }
  for (i in 151:200) {
    X[, i] <- -Z2 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 0
  }

  # 3. First White Noise Gap (201 to 300)
  for (i in 201:300) {
    X[, i] <- rnorm(n, mean = 0, sd = 1)
    beta_true[i] <- 0
  }

  # 4. Group 3: Inactive Trap 2 (301 to 360)
  for (i in 301:330) {
    X[, i] <- Z3 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 0
  }
  for (i in 331:360) {
    X[, i] <- -Z3 + rnorm(n, mean = 0, sd = sd_x)
    beta_true[i] <- 0
  }

  # 5. Second White Noise Gap (361 to p)
  for (i in 361:p) {
    X[, i] <- rnorm(n, mean = 0, sd = 1)
    beta_true[i] <- 0
  }

  # 6. Generate response y strictly via linear model
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X = X,
    y = y,
    beta_true = beta_true
  ))
}


# ------------------------------------------------------------------------------
# dgp_neg_traps_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP: active group + two inactive negative-correlation traps.
# Group 1 (1-100):   active,     +Z1/-Z1, beta = +3/-3.  s = 100.
# Group 2 (101-200): Trap 1,     +Z2/-Z2, beta = 0.
# Noise   (201-300): white noise.
# Group 3 (301-360): Trap 2,     +Z3/-Z3, beta = 0.
# Noise   (361-p):   white noise.
# sigma_y derived from target SNR.
#
# @param n     Number of observations.
# @param p     Number of predictors. Default 500.
# @param snr   Target signal-to-noise ratio.
# @param sd_x  Within-group noise sd. Default sqrt(0.01) => rho ~ 0.99.
# @param seed  Integer RNG seed. Default NULL.
#
# @return Named list: X, y, beta, n, p, s=100, snr, sigma_y, sd_x.

dgp_neg_traps_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  Z1 <- rnorm(n, 0, 1)
  Z2 <- rnorm(n, 0, 1)
  Z3 <- rnorm(n, 0, 1)

  X    <- matrix(0, nrow = n, ncol = p)
  beta <- rep(0, p)

  # Group 1: active (1-100)
  for (i in 1:50)   { X[, i] <- Z1 + rnorm(n, 0, sd_x);  beta[i] <- 3  }
  for (i in 51:100) { X[, i] <- -Z1 + rnorm(n, 0, sd_x); beta[i] <- -3 }

  # Group 2: Trap 1 (101-200)
  for (i in 101:150) X[, i] <- Z2 + rnorm(n, 0, sd_x)
  for (i in 151:200) X[, i] <- -Z2 + rnorm(n, 0, sd_x)

  # White noise (201-300)
  for (i in 201:300) X[, i] <- rnorm(n, 0, 1)

  # Group 3: Trap 2 (301-360)
  for (i in 301:330) X[, i] <- Z3 + rnorm(n, 0, sd_x)
  for (i in 331:360) X[, i] <- -Z3 + rnorm(n, 0, sd_x)

  # White noise (361-p)
  if (p >= 361L) for (i in 361:p) X[, i] <- rnorm(n, 0, 1)

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, 0, sigma_y)

  list(
    X       = X,
    y       = y,
    beta    = beta,
    n       = n,
    p       = p,
    s       = 100L,
    snr     = snr,
    sigma_y = sigma_y,
    sd_x    = sd_x
  )
}
