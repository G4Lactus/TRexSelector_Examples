# ==============================================================================
# dgp_nn_snr.R
#
# Data generating processes for the Nearest-Neighbours (NN) / MA(kappa).
#
# The covariance is banded via a causal MA(kappa) convolution:
#   X[i, j] = sum_{l=0}^{kappa} theta[l] * eta[i, j - l]
#
# with theta[l] = c * rho^l, normalised so that Var(X[i, j]) = 1.
#
# This yields Sigma[j, k] != 0 only for |j - k| <= kappa (banded structure).
#
# Support helpers (make_support_capped_spread, make_support_random) live in
# dgp_ar1_snr.R, which the demo file also sources.
# ==============================================================================

# ==============================================================================
# dgp_nn_snr
# ==============================================================================

#' Global MA(kappa) DGP with SNR-calibrated noise.
#'
#' Each column of X is the causal MA(kappa) convolution of a shared i.i.d.
#' N(0,1) innovation sheet eta:
#'
#'   X[i, j] = sum_{l=0}^{kappa} theta[l] * eta[i, j - l].
#'
#' In 1-based R indexing we avoid negative column indices by using a padded
#' innovation sheet eta of size n x (p + kappa). The same causal convolution is
#' then implemented as
#'
#'   X[, j] = eta[, j:(j + kappa)] %*% rev(theta),
#'
#' which is algebraically equivalent after reindexing.
#'
#' The coefficients are
#'
#'   theta[l] = c * rho^l,
#'   c = 1 / sqrt(sum_{l=0}^{kappa} rho^(2l)).
#'
#' The resulting covariance between columns j and k (d = |j - k|) is
#'
#'   Sigma[j, k] = sum_{l=0}^{kappa - d} theta[l] * theta[l + d],   d <= kappa,
#'   Sigma[j, k] = 0,                                                d >  kappa.
#'
#' The response is
#'
#'   y = X %*% beta + eps,   eps ~ N(0, noise_sigma^2 * I_n),
#'
#' where
#'
#'   noise_sigma = sqrt(sample_var(X %*% beta) / snr),
#'
#' using the usual sample variance with n - 1 denominator.
#'
#' @param n Number of observations.
#' @param p Number of predictors.
#' @param support 1-based integer vector of active predictor indices.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param kappa MA order (bandwidth). Default 3.
#' @param rho Decay parameter: theta[l] proportional to rho^l. Default 0.7.
#' @param snr Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param seed Integer seed for the RNG (governs both X and noise).
#'
#' @return A list with components:
#'   X, y, beta, true_support, n, p, s, kappa, rho, snr.
#'
#' @examples
#' support <- make_support_capped_spread(s = 5L, p = 500L, max_gap = 20L)
#' dat <- dgp_nn_snr(
#'   n = 150,
#'   p = 500,
#'   support = support,
#'   kappa = 3L,
#'   rho = 0.7,
#'   snr = 2.0,
#'   seed = 2026
#' )
dgp_nn_snr <- function(n, p, support,
                       amplitude = 3,
                       kappa = 3L,
                       rho = 0.7,
                       snr = 2.0,
                       seed = NULL) {

  stopifnot(kappa >= 1L)
  stopifnot(rho >= 0, rho < 1)
  stopifnot(snr > 0)
  stopifnot(all(support >= 1L), all(support <= p))

  if (!is.null(seed)) set.seed(seed)

  # Normalised MA coefficients:
  # theta[l] = c * rho^l, c = 1 / sqrt(sum_{l=0}^{kappa} rho^{2l})
  l_seq <- 0L:kappa
  theta <- rho ^ l_seq
  theta <- theta / sqrt(sum(theta^2))
  theta_rev <- rev(theta)

  # Innovation sheet: n x (p + kappa)
  eta <- matrix(stats::rnorm(n * (p + kappa)), nrow = n, ncol = p + kappa)

  # Build X via causal MA(kappa) convolution
  # Column j (1-based) = sum_{l=0}^kappa theta[l] * eta[, j + kappa - l]
  # Implemented as eta[, j:(j + kappa)] %*% rev(theta)
  X <- matrix(0.0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    X[, j] <- drop(eta[, j:(j + kappa)] %*% theta_rev)
  }

  # Coefficient vector (1-based support)
  beta <- numeric(p)
  beta[support] <- amplitude

  # SNR-calibrated response: noise_sigma = sqrt(var(signal) / snr)
  # var() uses the n - 1 denominator.
  signal <- drop(X %*% beta)
  var_sig <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(
    X = X,
    y = y,
    beta = beta,
    true_support = as.integer(support),
    n = n,
    p = p,
    s = length(support),
    kappa = kappa,
    rho = rho,
    snr = snr
  )
}

# ==============================================================================
# dgp_nn_block_snr
# ==============================================================================

#' Block-diagonal MA(kappa) DGP with SNR-calibrated noise.
#'
#' The p predictors are divided into n_blocks independent blocks, each of size
#' Q = p / n_blocks. Within block m the covariance is the same causal MA(kappa)
#' banded structure as in dgp_nn_snr(), restricted to Q predictors; between
#' blocks the correlation is exactly zero.
#'
#' This is the nearest-neighbour analogue of dgp_ar1_block_snr():
#' - dgp_ar1_block_snr => blocks with AR(1) / Toeplitz covariance
#' - dgp_nn_block_snr  => blocks with MA(kappa) banded covariance
#'
#' Each block draws its own independent innovation sheet, so block RNG draws
#' are sequential (block 1 first, then block 2, etc.). This differs from the
#' global dgp_nn_snr() where a single sheet covers all p + kappa columns.
#'
#' Active predictors are typically set to one per block via
#' make_support_bt_one_per_block() from dgp_bt_snr.R (deterministic: first
#' element of each block).
#'
#' @param n Number of observations.
#' @param p Number of predictors (must be divisible by n_blocks).
#' @param support 1-based integer vector of active predictor indices.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param n_blocks Number of MA(kappa) blocks (= number of active predictors
#'   when one active per block).
#' @param kappa MA order (bandwidth) within each block. Default 3.
#'   Must satisfy kappa < Q = p / n_blocks.
#' @param rho Decay parameter: theta[l] proportional to rho^l. Default 0.7.
#' @param snr Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param seed Integer seed for the RNG.
#'
#' @return A list with components:
#'   X, y, beta, true_support, n, p, s, n_blocks, kappa, rho, snr.
#'
#' @examples
#' support <- make_support_bt_one_per_block(p = 25L, n_blocks = 5L)
#' dat <- dgp_nn_block_snr(
#'   n = 150,
#'   p = 25,
#'   support = support,
#'   n_blocks = 5,
#'   kappa = 3L,
#'   rho = 0.7,
#'   snr = 2.0,
#'   seed = 2026
#' )
dgp_nn_block_snr <- function(n, p, support,
                             amplitude = 3,
                             n_blocks,
                             kappa = 3L,
                             rho = 0.7,
                             snr = 2.0,
                             seed = NULL) {

  stopifnot(p %% n_blocks == 0L)
  block_size <- p %/% n_blocks
  stopifnot(kappa < block_size,
            "kappa must be less than block size Q = p / n_blocks")
  stopifnot(kappa >= 1L)
  stopifnot(rho >= 0, rho < 1)
  stopifnot(snr > 0)
  stopifnot(all(support >= 1L), all(support <= p))

  if (!is.null(seed)) set.seed(seed)

  # Shared normalised MA coefficients (identical across all blocks)
  l_seq <- 0L:kappa
  theta <- rho ^ l_seq
  theta <- theta / sqrt(sum(theta^2))
  theta_rev <- rev(theta)

  # Build X block-by-block; each block has its own independent innovation sheet
  X <- matrix(0.0, nrow = n, ncol = p)

  for (m in seq_len(n_blocks)) {
    col_start <- (m - 1L) * block_size + 1L

    # Innovation sheet for this block: n x (block_size + kappa)
    eta_m <- matrix(stats::rnorm(n * (block_size + kappa)),
                    nrow = n, ncol = block_size + kappa)

    for (j in seq_len(block_size)) {
      X[, col_start + j - 1L] <- drop(eta_m[, j:(j + kappa)] %*% theta_rev)
    }
  }

  # Coefficient vector (1-based support)
  beta <- numeric(p)
  beta[support] <- amplitude

  # SNR-calibrated response
  signal <- drop(X %*% beta)
  var_sig <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(
    X = X,
    y = y,
    beta = beta,
    true_support = as.integer(support),
    n = n,
    p = p,
    s = length(support),
    n_blocks = n_blocks,
    kappa = kappa,
    rho = rho,
    snr = snr
  )
}
