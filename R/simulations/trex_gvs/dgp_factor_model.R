# ==============================================================================
# dgp_factor_model.R
# ==============================================================================
# Data-generating process for the overlapping latent factor model scenario.
# Sourced by demo_trex_gvs_08_4.R.
#
# Contents:
#   generate_factor_model_data(n, p, K, sigma_y)
#       Interactive demo DGP — overlapping latent factors with a small active subset.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_factor_model_data
# ------------------------------------------------------------------------------
# DGP Scenario: Realistic Factor Model.
# K = 4 latent factors driving the covariance matrix.
# Features have overlapping loadings (e.g., a gene belongs to two pathways).
# Only vars 1-20 are true actives (subset of Factor 1's range); rest are traps.

generate_factor_model_data <- function(n = 200, p = 500, K = 4, sigma_y = 15) {

  # 1. Generate Latent Factors (The unobserved drivers)
  F_factors <- matrix(rnorm(n * K), nrow = n, ncol = K)

  # 2. Define the Loading Matrix (Lambda)
  # Rows: Factors, Columns: p variables
  Lambda <- matrix(0, nrow = K, ncol = p)

  # Define overlapping "sectors" or "pathways"
  # Factor 1 drives vars 1 to 100
  Lambda[1, 1:100] <- runif(100, 0.5, 1.5)

  # Factor 2 drives vars 50 to 150 (Notice the overlap from 50-100!)
  Lambda[2, 50:150] <- runif(101, 0.5, 1.5)

  # Factor 3 drives vars 120 to 220
  Lambda[3, 120:220] <- runif(101, 0.5, 1.5)

  # Factor 4 drives vars 200 to 300
  Lambda[4, 200:300] <- runif(101, 0.5, 1.5)

  # Vars 301 to 500 have no factor loadings (pure independent noise)

  # 3. Generate Predictors X
  # X = F * Lambda + idiosyncratic noise
  idiosyncratic_noise <- matrix(rnorm(n * p, 0, sd = 0.5), nrow = n, ncol = p)
  X <- F_factors %*% Lambda + idiosyncratic_noise
  X <- scale(X) # Standardize for linear model assumptions

  # 4. Define True Support (beta)
  beta_true <- rep(0, p)

  # We define only a SUBSET of Factor 1's variables as active.
  # Variables 1 to 20 are true active variables.
  # Variables 21 to 100 are NULL traps, but they are highly correlated to 1-20
  # because they share Factor 1.
  # Variables 50 to 100 are even more complex traps because they also share Factor 2.

  active_indices <- 1:20
  beta_true[active_indices] <- 3.0

  # Randomize signs to ensure we test the absolute correlation logic
  signs <- sample(c(-1, 1), length(active_indices), replace = TRUE)
  beta_true[active_indices] <- beta_true[active_indices] * signs

  # 5. Generate response y strictly via linear model
  epsilon <- rnorm(n, mean = 0, sd = sigma_y)
  y <- as.numeric(X %*% beta_true + epsilon)
  y <- scale(y)[, 1]

  return(list(
    X = X,
    y = y,
    beta_true = beta_true
  ))
}


# TODO: dgp_factor_model_snr() — SNR-calibrated counterpart
