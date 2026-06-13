# ==============================================================================
# dgp_arma_blocks.R
# ==============================================================================
# Data-generating process for the ARMA-structured blocks scenario.
# Sourced by demo_trex_gvs_11.R.
#
# Contents:
#   generate_arma_blocks(n, p, ar_coef, ma_coefs, sigma_y)
#       Interactive demo DGP — 4 blocks governed by different ARMA processes.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_arma_blocks
# ------------------------------------------------------------------------------
# DGP Scenario: ARMA-Structured Blocks.
# Generates distinct blocks governed by different time-series processes.
# Block 1: AR(1)
# Block 2: MA(3)
# Block 3: ARMA(2,1)
# Block 4: Inactive Trap (AR(1))

generate_arma_blocks <- function(n = 200,
                                 p = 500,
                                 ar_coef = 0.8,
                                 ma_coefs = c(0.5, 0.3, 0.1),
                                 sigma_y = 15) {

  group_sizes <- c(20, 50, 80, 65)
  is_active <- c(TRUE, TRUE, TRUE, FALSE) # Block 4 (size 65) is the trap

  # Map blocks to specific ARMA models
  # Block 1: AR(1)
  model_1 <- list(ar = ar_coef)
  # Block 2: MA(3)
  model_2 <- list(ma = ma_coefs)
  # Block 3: ARMA(2,1)
  model_3 <- list(ar = c(ar_coef, ar_coef / 5), ma = ma_coefs[1])
  # Block 4 (Trap): AR(1) to see if EN swallows long-memory nulls
  model_4 <- list(ar = ar_coef)

  models <- list(model_1, model_2, model_3, model_4)

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
    # Place white noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rnorm(n * gaps[i])
      current_idx <- end_idx + 1
    }

    # Place the ARMA block
    if (i <= num_blocks) {
      b_id <- block_order[i]
      b_size <- group_sizes[b_id]
      end_idx <- current_idx + b_size - 1

      # Generate n independent realizations of the ARMA process of length b_size
      for (row in 1:n) {
        X[row, current_idx:end_idx] <- as.numeric(arima.sim(model = models[[b_id]], n = b_size))
      }

      # Standardize the block columns so variance is 1 (maintains EN scaling assumptions)
      X[, current_idx:end_idx] <- scale(X[, current_idx:end_idx])

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
  y <- scale(y)[, 1]

  return(list(
    X = X,
    y = y,
    beta_true = beta_true,
    block_indices = block_indices,
    block_order = block_order
  ))
}


# ==============================================================================
# dgp_arma_blocks_snr
# ==============================================================================
#' SNR-calibrated ARMA-structured blocks DGP.
#'
#' 4 blocks governed by different ARMA processes (column-standardized):
#'   Block 1 — AR(1);  Block 2 — MA(3);  Block 3 — ARMA(2,1);  Block 4 (trap) — AR(1).
#' 3 active blocks (sizes 20, 50, 80), 1 inactive trap (size 65).
#' Active set: s = 150.
#'
#' @param n       Number of observations.
#' @param p       Total predictors. Default 500.
#' @param snr     Signal-to-noise ratio.
#' @param ar_coef AR coefficient shared by blocks 1, 3, and the trap. Default 0.8.
#' @param ma_coefs MA coefficients for blocks 2 and 3. Default c(0.5, 0.3, 0.1).
#' @param seed    RNG seed. Default NULL.
#' @return list(X, y, beta, n, p, s, snr, sigma_y, ar_coef, block_indices, block_order)
dgp_arma_blocks_snr <- function(n, p = 500, snr, ar_coef = 0.8,
                                 ma_coefs = c(0.5, 0.3, 0.1),
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  group_sizes   <- c(20L, 50L, 80L, 65L)
  is_active     <- c(TRUE, TRUE, TRUE, FALSE)
  models        <- list(
    list(ar = ar_coef),
    list(ma = ma_coefs),
    list(ar = c(ar_coef, ar_coef / 5), ma = ma_coefs[1L]),
    list(ar = ar_coef)
  )
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
      for (row in seq_len(n)) {
        X[row, current_idx:end_idx] <-
          as.numeric(arima.sim(model = models[[b_id]], n = b_size))
      }
      # Standardize columns to variance 1
      X[, current_idx:end_idx] <- scale(X[, current_idx:end_idx])
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
    ar_coef       = ar_coef,
    block_indices = block_indices,
    block_order   = block_order
  )
}
