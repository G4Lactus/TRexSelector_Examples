# ==============================================================================
# dgp_chunked.R
# ==============================================================================
#
# Out-of-core data generating process for massive datasets.
# Generates a synthetic Gaussian design matrix X in column chunks and writes
# it directly to a binary file, bypassing R's memory limits.
#
# Mirrors the statistical properties of dgp_gauss_snr.R exactly.
#
# ==============================================================================

#' Chunked DGP for the i.i.d. Gaussian scenario with SNR-calibrated noise.
#'
#' Generates a massive design matrix X in memory-efficient column chunks,
#' writing it directly to a memory-mapped matrix using TRexSelector's native
#' `mmap_matrix` backend. Bypasses R's memory limits for large datasets.
#'
#' @param n            Number of observations.
#' @param p            Number of predictors.
#' @param support      1-based integer vector of active predictor indices.
#' @param x_path       File path where the binary matrix will be written.
#' @param chunk_cols   Number of columns to generate per RAM chunk (default 1000).
#' @param amplitude    Scalar signal coefficient applied to active predictors.
#' @param coefs        Optional numeric vector overriding amplitude.
#' @param snr          Target signal-to-noise ratio.
#' @param seed         Integer seed for the RNG.
#'
#' @return A list with components:
#'   \describe{
#'     \item{X}{Memory-mapped matrix object (TRexSelector mmap_matrix).}
#'     \item{y}{Response vector (length n, in memory).}
#'     \item{beta}{Coefficient vector (length p, in memory).}
#'     \item{true_support}{1-based integer index vector of active predictors.}
#'     \item{n, p, s, snr}{Dimension and parameter scalars.}
#'   }
dgp_chunked_gauss_snr <- function(n,
                                  p,
                                  support,
                                  x_path,
                                  chunk_cols = 1000L,
                                  amplitude  = 1,
                                  coefs      = NULL,
                                  snr        = 1.0,
                                  seed       = NULL) {

  stopifnot(is.integer(support) || is.numeric(support))
  stopifnot(all(support >= 1L), all(support <= p))
  stopifnot(snr > 0)

  if (!is.null(seed)) set.seed(seed)

  # Build coefficient vector in memory (p doubles is negligible for RAM)
  beta         <- numeric(p)
  active_coefs <- if (!is.null(coefs)) coefs else rep(amplitude, length(support))
  beta[support] <- active_coefs

  signal <- numeric(n)

  # Instantiate the memory-mapped matrix directly via TRexSelector
  X_mmap <- TRexSelector::mmap_matrix(x_path, n, p)

  # Generate and write X in chunks of `chunk_cols`
  start_col <- 1L
  while (start_col <= p) {
    end_col <- min(start_col + chunk_cols - 1L, p)
    ncols   <- end_col - start_col + 1L

    # Generate chunk
    X_chunk <- matrix(stats::rnorm(n * ncols), nrow = n, ncol = ncols)

    # Write chunk to the memory-mapped matrix
    X_mmap[, start_col:end_col] <- X_chunk

    # Accumulate signal chunk-by-chunk
    beta_chunk <- beta[start_col:end_col]
    signal <- signal + as.numeric(X_chunk %*% beta_chunk)

    start_col <- start_col + chunk_cols
  }

  # SNR-calibrated response (matches in-memory Bessel correction)
  var_sig     <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)

  # Note: This is drawn after all X random numbers, identical to the in-memory version.
  y <- signal + stats::rnorm(n, mean = 0, sd = noise_sigma)

  list(X = X_mmap, y = y, beta = beta, true_support = as.integer(support),
       n = n, p = p, s = length(support), snr = snr)
}
# ==============================================================================
