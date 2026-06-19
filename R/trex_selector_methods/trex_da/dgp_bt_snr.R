# ==============================================================================
# dgp_bt_snr.R
#
# Data generating process for the Binary-Tree (BT) / hierarchical-block
# scenario, using SNR-based noise calibration.
#
# Features:
#   - Noise is calibrated to a target SNR = Var(signal) / Var(noise).
#   - The true support is passed in explicitly (as 1-based R indices).
#
# Public functions:
# -------------------
# Factor-model block covariance:
#   dgp_bt_snr(n, p, support, amplitude, n_blocks, rho_within, rho_between,
#              snr, seed)
#
# AR(1) block covariance:
#   dgp_ar1_block_snr(n, p, support, amplitude, n_blocks, rho, snr, seed)
#
# M AR(1) blocks + 1 white-noise block:
#   dgp_ar1_block_white_snr(n, p_ar, p_white, support, amplitude, n_blocks,
#                           rho, snr, seed)
#
# ==============================================================================


# ==============================================================================
# dgp_bt_snr
# ==============================================================================

#' DGP for the binary-tree scenario with SNR-calibrated noise.
#'
#' Two-level hierarchical block model.
#' The p predictors are divided into n_blocks equal-sized blocks.
#' Each predictor is generated as:
#'
#'   X[i, j] = sqrt(rho_between)              * f_0[i]
#'           + sqrt(rho_within - rho_between) * f_m[i]   (m = block of j)
#'           + sqrt(1 - rho_within)           * eta[i, j]
#'
#' with f_0, f_m, eta all i.i.d. N(0, 1), yielding:
#'
#'   Sigma[j, k] = rho_within   (same block)
#'   Sigma[j, k] = rho_between  (different block)
#'
#' The response is:
#'
#'   y = X %*% beta + eps,   eps ~ N(0, noise_sigma^2 * I_n)
#'
#' where noise_sigma = sqrt( sample_var(X %*% beta) / snr ).
#'
#' @param n           Number of observations.
#' @param p           Number of predictors (must be divisible by n_blocks).
#' @param support     1-based integer vector of active predictor indices.
#' @param amplitude   Signal coefficient for each active predictor. Default 3.
#' @param n_blocks    Number of hierarchical blocks.
#' @param rho_within  Within-block correlation. Must satisfy rho_within > rho_between.
#' @param rho_between Between-block correlation. Must satisfy rho_between >= 0.
#' @param snr         Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param seed        Integer seed for the RNG (governs X and noise).
#'
#' @return A list with components:
#'   X, y, beta, true_support, n, p, s, n_blocks, rho_within, rho_between, snr.
#'
#' @examples
#' support <- make_support_capped_spread(s = 5L, p = 500L, max_gap = 20L)
#' dat <- dgp_bt_snr(
#'   n           = 150,
#'   p           = 500,
#'   support     = support,
#'   amplitude   = 3,
#'   n_blocks    = 10,
#'   rho_within  = 0.8,
#'   rho_between = 0.1,
#'   snr         = 2.0,
#'   seed        = 2026
#' )
dgp_bt_snr <- function(n,
                       p,
                       support,
                       amplitude   = 3,
                       n_blocks    = 10L,
                       rho_within  = 0.8,
                       rho_between = 0.1,
                       snr         = 2.0,
                       seed        = NULL) {

  stopifnot(p %% n_blocks == 0)
  stopifnot(rho_between >= 0, rho_between <= rho_within, rho_within < 1)
  stopifnot(snr > 0)
  stopifnot(all(support >= 1L), all(support <= p))

  if (!is.null(seed)) set.seed(seed)

  block_size   <- p %/% n_blocks
  block_labels <- rep(seq_len(n_blocks), each = block_size)  # length p

  a  <- sqrt(rho_between)
  b  <- sqrt(rho_within - rho_between)
  c <- sqrt(1.0 - rho_within)

  # Factor draws:
  #   row i: 1 global draw (f0), n_blocks block draws, p idiosyncratic draws.
  f0      <- stats::rnorm(n)
  F_block <- matrix(stats::rnorm(n * n_blocks), nrow = n, ncol = n_blocks)
  eta     <- matrix(stats::rnorm(n * p),        nrow = n, ncol = p)

  X <- a  * matrix(f0, nrow = n, ncol = p) +
    b  * F_block[, block_labels, drop = FALSE] +
    c * eta

  # Coefficient vector
  beta          <- numeric(p)
  beta[support] <- amplitude

  # SNR-calibrated response: noise_sigma = sqrt(var(signal) / snr)
  signal      <- drop(X[, support, drop = FALSE] %*% beta[support])
  var_sig     <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y           <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(
    X            = X,
    y            = y,
    beta         = beta,
    true_support = as.integer(support),
    n            = n,
    p            = p,
    s            = length(support),
    n_blocks     = n_blocks,
    rho_within   = rho_within,
    rho_between  = rho_between,
    snr          = snr
  )
}



# ==============================================================================
# dgp_ar1_block_snr
# ==============================================================================

#' DGP for the paper's block-diagonal AR(1)/Toeplitz scenario.
#'
#' The p predictors are divided into n_blocks equal-sized blocks.
#' Within each block the predictors follow an AR(1) process with parameter rho:
#'
#'   X[i, block_start]     ~ N(0, 1)
#'   X[i, block_start + j] = rho * X[i, block_start + j - 1]
#'                         + sqrt(1 - rho^2) * eps[i, j],  eps ~ N(0,1)
#'
#' yielding the Toeplitz covariance Sigma_m[j, k] = rho^|j-k| within each
#' block and zero correlation between blocks (block-diagonal Sigma).
#'
#'
#' @param n         Number of observations.
#' @param p         Number of predictors (must be divisible by n_blocks).
#' @param support   1-based integer vector of active predictor indices.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param n_blocks  Number of AR(1) blocks.
#' @param rho       AR(1) parameter (within-block). Must satisfy |rho| < 1.
#' @param snr       Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param seed      Integer seed for the RNG.
#'
#' @return A list with components:
#'   X, y, beta, true_support, n, p, s, n_blocks, rho, snr.
dgp_ar1_block_snr <- function(n,
                              p,
                              support,
                              amplitude = 3,
                              n_blocks,
                              rho,
                              snr  = 2.0,
                              seed = NULL) {

  stopifnot(p %% n_blocks == 0)
  stopifnot(abs(rho) < 1)
  stopifnot(snr > 0)
  stopifnot(all(support >= 1L), all(support <= p))

  if (!is.null(seed)) set.seed(seed)

  block_size <- p %/% n_blocks
  innov_sd   <- sqrt(1.0 - rho^2)

  # Build X block by block using the AR(1) column recursion.
  X <- matrix(0.0, nrow = n, ncol = p)
  for (m in seq_len(n_blocks)) {
    col_start <- (m - 1L) * block_size + 1L
    col_end   <- m * block_size

    X[, col_start] <- stats::rnorm(n)
    for (j in seq(col_start + 1L, col_end)) {
      X[, j] <- rho * X[, j - 1L] + innov_sd * stats::rnorm(n)
    }
  }

  # Coefficient vector
  beta          <- numeric(p)
  beta[support] <- amplitude

  # SNR-calibrated response: noise_sigma = sqrt(var(signal) / snr)
  signal      <- drop(X[, support, drop = FALSE] %*% beta[support])
  var_sig     <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y           <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(
    X            = X,
    y            = y,
    beta         = beta,
    true_support = as.integer(support),
    n            = n,
    p            = p,
    s            = length(support),
    n_blocks     = n_blocks,
    rho          = rho,
    snr          = snr
  )
}



# ==============================================================================
# dgp_ar1_block_white_snr
# ==============================================================================

#' DGP: M times AR(1) blocks + 1 white-noise block; actives only in the AR blocks.
#'
#' The first p_ar predictors are divided into n_blocks equal-sized AR(1) blocks.
#' The last p_white predictors form one independent white-noise block:
#' each column is i.i.d. N(0, 1), uncorrelated with all other columns.
#'
#' Active predictors (support) must all lie within the AR part (indices 1 ..
#' p_ar).
#' The selector operates on the full design matrix of total width:
#' p = p_ar + p_white.
#'
#' @param n         Number of observations.
#' @param p_ar      Number of AR(1) predictors (must be divisible by n_blocks).
#' @param p_white   Number of white-noise predictors (>= 1).
#' @param support   1-based integer vector of active indices; all <= p_ar.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param n_blocks  Number of AR(1) blocks.
#' @param rho       AR(1) parameter (within AR blocks). Must satisfy |rho| < 1.
#' @param snr       Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param seed      Integer seed for the RNG.
#'
#' @return A list with components:
#'   X, y, beta, true_support, n, p (= p_ar + p_white), p_ar, p_white,
#'   s, n_blocks, rho, snr.
dgp_ar1_block_white_snr <- function(n,
                                    p_ar,
                                    p_white,
                                    support,
                                    amplitude = 3,
                                    n_blocks,
                                    rho,
                                    snr  = 2.0,
                                    seed = NULL) {

  stopifnot(p_ar %% n_blocks == 0)
  stopifnot(p_white >= 1L)
  stopifnot(abs(rho) < 1)
  stopifnot(snr > 0)
  stopifnot(all(support >= 1L), all(support <= p_ar))

  if (!is.null(seed)) set.seed(seed)

  p          <- p_ar + p_white
  block_size <- p_ar %/% n_blocks
  innov_sd   <- sqrt(1.0 - rho^2)

  # Build AR part block by block.
  X_ar <- matrix(0.0, nrow = n, ncol = p_ar)
  for (m in seq_len(n_blocks)) {
    col_start <- (m - 1L) * block_size + 1L
    col_end   <- m * block_size

    X_ar[, col_start] <- stats::rnorm(n)
    for (j in seq(col_start + 1L, col_end)) {
      X_ar[, j] <- rho * X_ar[, j - 1L] + innov_sd * stats::rnorm(n)
    }
  }

  # White-noise block: i.i.d. N(0, 1) — independent of everything.
  X_white <- matrix(stats::rnorm(n * p_white), nrow = n, ncol = p_white)

  X <- cbind(X_ar, X_white)

  # Coefficient vector (1-based support, restricted to AR part)
  beta          <- numeric(p)
  beta[support] <- amplitude

  # SNR-calibrated response
  signal      <- drop(X %*% beta)
  var_sig     <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y           <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(
    X            = X,
    y            = y,
    beta         = beta,
    true_support = as.integer(support),
    n            = n,
    p            = p,
    p_ar         = p_ar,
    p_white      = p_white,
    s            = length(support),
    n_blocks     = n_blocks,
    rho          = rho,
    snr          = snr
  )
}
