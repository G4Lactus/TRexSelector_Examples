# ==============================================================================
# dgp_hastie.R
# ==============================================================================
# Data-generating processes for the Hastie (Zou & Hastie 2005) GVS scenario.
# Sourced by demo_trex_gvs_01.R and simulation_trex_gvs_01.R.
#
# Contents:
#   generate_hastie_scaled_data(n, p, sigma_y, sd_x)
#       Interactive demo DGP — direct noise parameter, used for single-run plots.
#   dgp_hastie_snr(n, p, snr, sd_x, seed)
#       Simulation-ready DGP — SNR-calibrated noise, uniform return contract.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_hastie_scaled_data
# ------------------------------------------------------------------------------
# Data generation model from Zou and Hastie (2005) - Example 4
# "Regularization and variable selection via the elastic net"
# 3 groups of 50 equi-correlated active predictors (indices 1-150); rest noise.

generate_hastie_scaled_data <- function(n = 100, p = 500, sigma_y = 15,
                                        sd_x = sqrt(0.01)) {

  if (p < 150) {
    stop("Error: p must be at least 150 to accommodate the 3 groups.")
  }

  # 1. Define the true beta coefficient vector (The Support Vector)
  # Hastie's rule: ALL variables in the correlated groups are true actives.
  beta_true <- rep(0, p)
  beta_true[1:150] <- 3

  # Initialize the predictor matrix X
  X <- matrix(0, nrow = n, ncol = p)

  # 2. Generate latent driving factors for the 3 groups
  Z1 <- rnorm(n, mean = 0, sd = 1)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  Z3 <- rnorm(n, mean = 0, sd = 1)

  # 3. Generate the 150 grouped predictors
  # Group 1 (Variables 1 to 50)
  for (i in 1:50) {
    X[, i] <- Z1 + rnorm(n, mean = 0, sd = sd_x)
  }
  # Group 2 (Variables 51 to 100)
  for (i in 51:100) {
    X[, i] <- Z2 + rnorm(n, mean = 0, sd = sd_x)
  }
  # Group 3 (Variables 101 to 150)
  for (i in 101:150) {
    X[, i] <- Z3 + rnorm(n, mean = 0, sd = sd_x)
  }

  # 4. Generate the remainder as independent white noise
  if (p > 150) {
    for (i in 151:p) {
      X[, i] <- rnorm(n, mean = 0, sd = 1)
    }
  }

  # 5. Generate the response variable y
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)

  return(list(
    X         = X,
    y         = y,
    beta_true = beta_true
  ))
}


# ------------------------------------------------------------------------------
# dgp_hastie_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP for the Hastie GVS scenario with SNR-calibrated noise.
#
# The within-group correlation is rho = 1 / (1 + sd_x^2).
# To sweep rho, convert: sd_x = sqrt((1 - rho) / rho).
#
# The noise standard deviation sigma_y is derived from the target SNR:
#
#   signal   = X %*% beta
#   sigma_y  = sqrt( var(signal) / snr )      (n-1 denominator)
#
# @param n     Number of observations.
# @param p     Number of predictors (>= 150). Default 500.
# @param snr   Target signal-to-noise ratio Var(signal) / sigma_y^2.
# @param sd_x  Within-group noise sd. Controls rho = 1/(1+sd_x^2).
#              Default sqrt(0.01) => rho ~ 0.99.
# @param seed  Integer RNG seed. Default NULL.
#
# @return Named list:
#   X        n x p predictor matrix
#   y        length-n response
#   beta     length-p true coefficient vector (field consistent with .run_mc)
#   n, p, s  dimensions
#   snr      target SNR
#   sigma_y  realised noise standard deviation
#   sd_x     within-group noise sd used

dgp_hastie_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01), seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (p < 150) {
    stop("'p' must be at least 150 to accommodate the 3 groups.")
  }

  # True coefficient vector
  beta <- rep(0, p)
  beta[1:150] <- 3

  # Generate predictors via three latent factors
  X <- matrix(0, nrow = n, ncol = p)

  Z1 <- rnorm(n, mean = 0, sd = 1)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  Z3 <- rnorm(n, mean = 0, sd = 1)

  for (i in 1:50)    X[, i] <- Z1 + rnorm(n, mean = 0, sd = sd_x)
  for (i in 51:100)  X[, i] <- Z2 + rnorm(n, mean = 0, sd = sd_x)
  for (i in 101:150) X[, i] <- Z3 + rnorm(n, mean = 0, sd = sd_x)
  for (i in 151:p)   X[, i] <- rnorm(n, mean = 0, sd = 1)

  # SNR-calibrated noise: sigma_y = sqrt(var(signal) / snr)
  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)    # var() uses n-1 denominator

  # Generate response
  y <- signal + rnorm(n, mean = 0, sd = sigma_y)

  list(
    X       = X,
    y       = y,
    beta    = beta,
    n       = n,
    p       = p,
    s       = 150L,
    snr     = snr,
    sigma_y = sigma_y,
    sd_x    = sd_x
  )
}
