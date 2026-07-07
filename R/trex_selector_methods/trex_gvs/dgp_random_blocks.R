# ==============================================================================
# dgp_random_blocks.R
# ==============================================================================
# Data-generating process for the randomly ordered, unequal separated blocks scenario.
# Not currently sourced by any demo; kept for future scenarios.
#
# Contents:
#   generate_random_blocks_data(n, p, sigma_y)
#       Interactive demo DGP — 3 active blocks, shuffled with random noise gaps.
# ==============================================================================


# ------------------------------------------------------------------------------
# generate_random_blocks_data
# ------------------------------------------------------------------------------
# DGP Scenario: Randomly Ordered, Unequal Separated Blocks.
# Three blocks (sizes 20, 50, 80).
# Their order is shuffled randomly, and they are dropped into the p-dimensional
# space separated by random blocks of white noise.

generate_random_blocks_data <- function(n = 200, p = 500, sigma_y = 15) {

  group_sizes <- c(20, 50, 80)
  num_groups <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  if (p < total_grouped + num_groups - 1) {
    stop("Error: p is too small to hold groups and separating noise.")
  }

  # 1. Shuffle the order of the groups
  group_order <- sample(1:num_groups)
  shuffled_sizes <- group_sizes[group_order]

  # 2. Distribute the noise into (num_groups + 1) gaps
  num_noise <- p - total_grouped
  rem_noise <- num_noise - (num_groups - 1) # Reserve minimum 1 var for internal gaps

  # Randomly allocate remaining noise across all gaps
  gaps <- as.numeric(rmultinom(1, size = rem_noise, prob = rep(1, num_groups + 1)))

  # Add the reserved 1 to the internal gaps to guarantee separation
  for(i in 2:num_groups) {
    gaps[i] <- gaps[i] + 1
  }

  # Initialize
  X <- matrix(0, nrow = n, ncol = p)
  beta_true <- rep(0, p)
  sd_x <- sqrt(0.01)

  # Latent driving factors for the 3 original groups
  latent_factors <- list(
    rnorm(n, mean = 0, sd = 1),
    rnorm(n, mean = 0, sd = 1),
    rnorm(n, mean = 0, sd = 1)
  )

  # 3. Build the matrix sequentially
  current_idx <- 1
  group_indices <- list()

  for (i in 1:(num_groups + 1)) {
    # Place noise gap
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1
      X[, current_idx:end_idx] <- rnorm(n * gaps[i], mean = 0, sd = 1)
      current_idx <- end_idx + 1
    }

    # Place the block (if not the last outer gap)
    if (i <= num_groups) {
      block_id <- group_order[i]
      block_size <- shuffled_sizes[i]
      end_idx <- current_idx + block_size - 1

      Z <- latent_factors[[block_id]]
      for (c in current_idx:end_idx) {
        X[, c] <- Z + rnorm(n, mean = 0, sd = sd_x)
      }

      beta_true[current_idx:end_idx] <- 3
      group_indices[[block_id]] <- current_idx:end_idx

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
    group_indices = group_indices,
    group_order = group_order
  ))
}


# ------------------------------------------------------------------------------
# dgp_random_blocks_snr
# ------------------------------------------------------------------------------
# Simulation-ready DGP: 3 active blocks (20, 50, 80), randomly ordered and
# placed in the p-dimensional space with random noise gaps between them.
# Active set: s = 150. Within-block correlation rho = 1 / (1 + sd_x^2).
# sigma_y derived from target SNR.
#
# @param n     Number of observations.
# @param p     Number of predictors. Default 500.
# @param snr   Target signal-to-noise ratio.
# @param sd_x  Within-block noise sd. Default sqrt(0.01) => rho ~ 0.99.
# @param seed  Integer RNG seed. Default NULL.
#
# @return Named list: X, y, beta, n, p, s=150, snr, sigma_y, sd_x,
#                     group_sizes, group_order, group_indices.

dgp_random_blocks_snr <- function(n, p = 500, snr, sd_x = sqrt(0.01),
                                  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  group_sizes   <- c(20L, 50L, 80L)
  num_groups    <- length(group_sizes)
  total_grouped <- sum(group_sizes)

  group_order    <- sample(seq_len(num_groups))
  shuffled_sizes <- group_sizes[group_order]
  num_noise      <- p - total_grouped
  rem_noise      <- num_noise - (num_groups - 1L)
  gaps           <- as.numeric(rmultinom(1, size = rem_noise,
                                         prob = rep(1, num_groups + 1)))
  for (i in 2:num_groups) gaps[i] <- gaps[i] + 1

  X             <- matrix(0, nrow = n, ncol = p)
  beta          <- rep(0, p)
  latent        <- lapply(seq_len(num_groups), function(k) rnorm(n, 0, 1))
  current_idx   <- 1L
  group_indices <- vector("list", num_groups)

  for (i in 1:(num_groups + 1L)) {
    if (gaps[i] > 0) {
      end_idx <- current_idx + gaps[i] - 1L
      X[, current_idx:end_idx] <- matrix(rnorm(n * gaps[i], 0, 1), n, gaps[i])
      current_idx <- end_idx + 1L
    }
    if (i <= num_groups) {
      block_id   <- group_order[i]
      block_size <- shuffled_sizes[i]
      end_idx    <- current_idx + block_size - 1L
      Z          <- latent[[block_id]]
      for (c in current_idx:end_idx) X[, c] <- Z + rnorm(n, 0, sd_x)
      beta[current_idx:end_idx]  <- 3
      group_indices[[block_id]]  <- current_idx:end_idx
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
    s             = as.integer(total_grouped),
    snr           = snr,
    sigma_y       = sigma_y,
    sd_x          = sd_x,
    group_sizes   = group_sizes,
    group_order   = group_order,
    group_indices = group_indices
  )
}
