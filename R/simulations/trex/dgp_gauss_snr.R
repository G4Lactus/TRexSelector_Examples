# ==============================================================================
# dgp_gauss_snr.R
# ==============================================================================
#
# Data generating process for the standard i.i.d. Gaussian scenario, using
# SNR-based noise calibration.
#
# Mirrors C++ datagen::SyntheticData with predictor_policy::Normal() and
# noisegen::noise_policy::Normal().
#
# Public functions:
#   dgp_gauss_snr(n, p, support, amplitude, coefs, snr, seed)
#
# ==============================================================================


# ==============================================================================
# dgp_gauss_snr
# ==============================================================================

#' DGP for the i.i.d. Gaussian scenario with SNR-calibrated noise.
#'
#' Predictors are drawn i.i.d. from N(0, 1):
#'
#'   X[i, j]  i.i.d. N(0, 1)
#'
#' The response is:
#'
#'   y = X %*% beta + eps,   eps ~ N(0, noise_sigma^2 * I_n)
#'
#' where beta has entries equal to `amplitude` (or `coefs`) at the positions
#' given by `support` (1-based) and zero elsewhere.
#'
#' The noise standard deviation is chosen according to the target SNR:
#'
#'   noise_sigma = sqrt( var(X %*% beta) / snr )
#'
#' with var() using the n-1 denominator (Bessel correction), matching
#' C++ noisegen::calculate_noise_params().
#'
#' @param n         Number of observations.
#' @param p         Number of predictors.
#' @param support   1-based integer vector of active predictor indices.
#' @param amplitude Scalar signal coefficient applied to all active predictors.
#'                  Ignored when `coefs` is supplied. Default 1.
#' @param coefs     Optional numeric vector of length `length(support)` with
#'                  per-predictor active coefficients. Overrides `amplitude`.
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
#' dat <- dgp_gauss_snr(
#'   n       = 300,
#'   p       = 1000,
#'   support = c(1L, 21L, 41L, 61L, 81L, 101L, 121L, 141L, 161L, 181L),
#'   snr     = 1.0,
#'   seed    = 42L
#' )
dgp_gauss_snr <- function(n,
                           p,
                           support,
                           amplitude = 1,
                           coefs     = NULL,
                           snr       = 1.0,
                           seed      = NULL) {

  stopifnot(is.integer(support) || is.numeric(support))
  stopifnot(all(support >= 1L), all(support <= p))
  stopifnot(snr > 0)

  if (!is.null(seed)) set.seed(seed)

  # Generate X: n x p matrix of i.i.d. N(0, 1)
  # Drawn column-major (rnorm(n * p) fills column by column in R).
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)

  # Build coefficient vector
  beta         <- numeric(p)
  active_coefs <- if (!is.null(coefs)) coefs else rep(amplitude, length(support))
  beta[support] <- active_coefs

  # SNR-calibrated response
  signal      <- drop(X %*% beta)
  var_sig     <- stats::var(signal)   # uses n-1 denominator (Bessel correction)
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
    snr          = snr
  )
}
# ==============================================================================
