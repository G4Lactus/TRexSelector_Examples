# ==============================================================================
# dgp_ar1_snr.R
#
# Data generating process for the AR(1) Toeplitz scenario, using SNR-based
# noise calibration:
#
# - SNR = Var(signal) / Var(noise).
#
# Public functions:
# - dgp_ar1_snr(n, p, support, amplitude, rho, snr, seed)
#
# ==============================================================================

# ==============================================================================
# dgp_ar1_snr
# ==============================================================================

#' DGP for the AR(1) scenario with SNR-calibrated noise.
#'
#' Predictor rows are generated column-by-column via the stationary AR(1)
#' recursion:
#'
#'   X[i, 1]   = eta[i, 1]
#'   X[i, j]   = rho * X[i, j-1] + sqrt(1 - rho^2) * eta[i, j]
#'
#' with eta[i, j] i.i.d. N(0, 1).
#' This yields the Toeplitz covariance Sigma[j, k] = rho^|j - k|.
#'
#' The response is:
#'
#'   y = X %*% beta + eps,   eps ~ N(0, noise_sigma^2 * I_n)
#'
#' where beta has entries equal to `amplitude` at the positions given by
#' `support` (1-based) and zero elsewhere.
#' The noise standard deviation is chosen according to the target SNR:
#'
#'   noise_sigma = sqrt( sample_var(X %*% beta) / snr )
#'
#' with sample_var() using the n-1 denominator, matching C++ make_snr_response().
#'
#' @param n         Number of observations.
#' @param p         Number of predictors.
#' @param support   1-based integer vector of active predictor indices.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param rho       AR(1) correlation parameter rho in (-1, 1).
#' @param snr       Target signal-to-noise ratio  Var(signal) / Var(noise).
#' @param seed      Integer seed for the RNG (governs both X and noise).
#'
#' @return A list with components:
#'   \describe{
#'     \item{X}{Predictor matrix (n x p).}
#'     \item{y}{Response vector (length n).}
#'     \item{beta}{Coefficient vector (length p).}
#'     \item{true_support}{1-based integer index vector of active predictors.}
#'     \item{n, p, s, snr}{Dimension and parameter scalars.}
#'   }
#'
#' @examples
#' dat <- dgp_ar1_snr(
#'   n         = 150,
#'   p         = 500,
#'   support   = make_support_random(5, 500, seed = 2026),
#'   amplitude = 3,
#'   rho       = 0.7,
#'   snr       = 2.0,
#'   seed      = 2026
#' )
dgp_ar1_snr <- function(n,
                        p,
                        support,
                        amplitude = 3,
                        rho       = 0.7,
                        snr       = 2.0,
                        seed      = NULL) {

  stopifnot(is.integer(support) || is.numeric(support))
  stopifnot(all(support >= 1L), all(support <= p))
  stopifnot(abs(rho) < 1)
  stopifnot(snr > 0)

  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------------------------
  # Generate X: AR(1) recursion row-wise
  # eta[i, j] are drawn in column order: col 1, col 2, ..., col p
  # (each rnorm(n) call advances the RNG by n draws)
  # ------------------------------------------------------------------
  scale_ar <- sqrt(1.0 - rho^2)

  X      <- matrix(0.0, nrow = n, ncol = p)
  X[, 1] <- stats::rnorm(n)
  for (j in seq(2L, p)) {
    X[, j] <- rho * X[, j - 1L] + scale_ar * stats::rnorm(n)
  }

  # ------------------------------------------------------------------
  # Coefficient vector (1-based support)
  # ------------------------------------------------------------------
  beta          <- numeric(p)
  beta[support] <- amplitude

  # ------------------------------------------------------------------
  # SNR-calibrated response:
  # noise_sigma = sqrt(var(signal) / snr)
  # var() in R uses n-1 denominator (Bessel correction factor)
  # ------------------------------------------------------------------

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
    s            = length(support),
    snr          = snr,
    rho          = rho
  )
}
