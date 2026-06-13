# ==============================================================================
# dgp_ht_snr.R
#
# Data generating processes for heavy-tailed (Student-t) predictor matrices
# with block-diagonal AR(1) covariance structure, using SNR-based noise
# calibration consistent with dgp_bt_snr.R and dgp_ar1_snr.R.
#
# Each predictor block is generated as:
#   1. A Gaussian AR(1) block  Z ~ N(0, Sigma_AR1)
#   2. Divided element-wise by the square root of an independent chi-squared
#      scale factor to obtain a multivariate Student-t distribution with nu df.
#
# Two variants:
#   dgp_ht_block_snr()       — All p columns form M heavy-tailed AR(1) blocks.
#   dgp_ht_block_white_snr() — M heavy-tailed AR(1) blocks + p_white
#                              heavy-tailed i.i.d. (identity) columns.
#                              Allows p_total to remain fixed while M and Q vary.
#
# Both support two noise scenarios:
#   heavy_noise = FALSE  — Gaussian noise
#   heavy_noise = TRUE   — Student-t noise with the same nu
# ==============================================================================


# ==============================================================================
# dgp_ht_block_snr
# ==============================================================================

#' DGP for heavy-tailed AR(1)-block data with SNR-calibrated noise.
#'
#' All p predictors form M independent AR(1) blocks of size p/M.  Each block
#' is drawn as a multivariate Student-t via the Gaussian-scale-mixture
#' representation: Z ~ N(0, Sigma_AR1), then X = Z / sqrt(U/nu) where
#' U ~ chi2(nu) is shared across all columns within each observation.
#'
#' @param n           Number of observations.
#' @param p           Total number of predictors. Must be divisible by n_blocks.
#' @param support     1-based integer vector of active predictor indices.
#' @param amplitude   Signal coefficient for each active predictor. Default 1.
#' @param n_blocks    Number of AR(1) blocks (= number of active predictors).
#' @param rho         AR(1) within-block correlation parameter.
#' @param snr         Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param nu          Degrees of freedom for the Student-t distribution. Default 3.
#' @param heavy_noise If TRUE, noise is also Student-t(nu); otherwise Gaussian.
#' @param seed        Integer RNG seed. NULL means do not set seed.
#'
#' @return A list with components: X, y, beta, true_support, n, p, rho, snr,
#'   nu, heavy_noise.
dgp_ht_block_snr <- function(n,
                             p,
                             support,
                             amplitude   = 1.0,
                             n_blocks    = 5L,
                             rho         = 0.8,
                             snr         = 2.0,
                             nu          = 3.0,
                             heavy_noise = FALSE,
                             seed        = NULL) {

  if (!is.null(seed)) set.seed(seed)

  block_size <- p %/% n_blocks
  innov_sd   <- sqrt(1.0 - rho^2)

  # Step 1: Generate underlying Gaussian AR(1) blocks Z ~ N(0, Sigma_AR1)
  Z <- matrix(0.0, nrow = n, ncol = p)
  for (m in seq_len(n_blocks)) {
    col_start <- (m - 1L) * block_size + 1L
    col_end   <- m * block_size

    Z[, col_start] <- stats::rnorm(n)
    for (j in seq(col_start + 1L, col_end)) {
      Z[, j] <- rho * Z[, j - 1L] + innov_sd * stats::rnorm(n)
    }
  }

  # Step 2: Convert to multivariate Student-t via Gaussian scale mixture.
  # One chi-squared scale draw per observation, shared across all columns.
  U            <- stats::rchisq(n, df = nu)
  scale_factor <- sqrt(U / nu)
  X            <- Z / scale_factor

  # Step 3: Response and SNR-calibrated noise
  beta    <- numeric(p)
  beta[support] <- amplitude

  signal      <- drop(X %*% beta)
  noise_sigma <- sqrt(stats::var(signal) / snr)

  if (heavy_noise) {
    raw_noise <- stats::rt(n, df = nu)
    raw_var   <- nu / (nu - 2)
    eps <- raw_noise * (noise_sigma / sqrt(raw_var))
  } else {
    eps <- stats::rnorm(n, mean = 0, sd = noise_sigma)
  }

  y <- signal + eps

  list(
    X           = X,
    y           = y,
    beta        = beta,
    true_support = as.integer(support),
    n           = n,
    p           = p,
    rho         = rho,
    snr         = snr,
    nu          = nu,
    heavy_noise = heavy_noise
  )
}


# ==============================================================================
# dgp_ht_block_white_snr
# ==============================================================================

#' DGP for heavy-tailed AR(1)-block + identity-block data with SNR-calibrated noise.
#'
#' The p_total = p_ar + p_white predictors consist of:
#'   - p_ar columns: M independent AR(1) blocks of size p_ar/M, each drawn
#'     as heavy-tailed (Student-t) via the Gaussian scale mixture.
#'   - p_white columns: a single heavy-tailed i.i.d. block (identity covariance).
#'
#' Active predictors lie entirely within the AR blocks (one per block).
#' The white-noise columns allow p_total to stay fixed while M and Q vary.
#'
#' @param n           Number of observations.
#' @param p_ar        Number of AR(1) predictors. Must be divisible by n_blocks.
#' @param p_white     Number of white-noise predictors.
#' @param support     1-based integer vector of active predictor indices
#'                    (into the full p_ar + p_white column space).
#' @param amplitude   Signal coefficient for each active predictor. Default 1.
#' @param n_blocks    Number of AR(1) blocks.
#' @param rho         AR(1) within-block correlation parameter.
#' @param snr         Target signal-to-noise ratio Var(signal) / Var(noise).
#' @param nu          Degrees of freedom for the Student-t distribution. Default 3.
#' @param heavy_noise If TRUE, noise is also Student-t(nu); otherwise Gaussian.
#' @param seed        Integer RNG seed. NULL means do not set seed.
#'
#' @return A list with components: X, y, beta, true_support, n, p (= p_ar + p_white),
#'   rho, snr, nu, heavy_noise.
dgp_ht_block_white_snr <- function(n,
                                   p_ar,
                                   p_white,
                                   support,
                                   amplitude   = 1.0,
                                   n_blocks    = 5L,
                                   rho         = 0.8,
                                   snr         = 2.0,
                                   nu          = 3.0,
                                   heavy_noise = FALSE,
                                   seed        = NULL) {
  if (!is.null(seed)) set.seed(seed)

  p          <- p_ar + p_white
  block_size <- p_ar %/% n_blocks
  innov_sd   <- sqrt(1.0 - rho^2)

  X <- matrix(0.0, nrow = n, ncol = p)

  # AR(1) blocks with heavy-tailed scaling (one chi-sq draw per observation per block)
  for (m in seq_len(n_blocks)) {
    col_start <- (m - 1L) * block_size + 1L
    col_end   <- m * block_size

    Z_block        <- matrix(0.0, nrow = n, ncol = block_size)
    Z_block[, 1L]  <- stats::rnorm(n)
    for (j in seq(2L, block_size)) {
      Z_block[, j] <- rho * Z_block[, j - 1L] + innov_sd * stats::rnorm(n)
    }

    U_block            <- stats::rchisq(n, df = nu)
    scale_factor_block <- sqrt(U_block / nu)
    X[, col_start:col_end] <- Z_block / scale_factor_block
  }

  # White-noise block with heavy-tailed scaling (one chi-sq draw per observation)
  if (p_white > 0L) {
    Z_white            <- matrix(stats::rnorm(n * p_white), nrow = n, ncol = p_white)
    U_white            <- stats::rchisq(n, df = nu)
    scale_factor_white <- sqrt(U_white / nu)
    X[, (p_ar + 1L):p] <- Z_white / scale_factor_white
  }

  beta <- numeric(p)
  beta[support] <- amplitude

  signal      <- drop(X %*% beta)
  noise_sigma <- sqrt(stats::var(signal) / snr)

  if (heavy_noise) {
    raw_noise <- stats::rt(n, df = nu)
    raw_var   <- nu / (nu - 2)
    eps <- raw_noise * (noise_sigma / sqrt(raw_var))
  } else {
    eps <- stats::rnorm(n, mean = 0, sd = noise_sigma)
  }

  y <- signal + eps

  list(
    X           = X,
    y           = y,
    beta        = beta,
    true_support = as.integer(support),
    n           = n,
    p           = p,
    rho         = rho,
    snr         = snr,
    nu          = nu,
    heavy_noise = heavy_noise
  )
}
