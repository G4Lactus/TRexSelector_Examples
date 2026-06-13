# ==============================================================================
# dgp_hapgen_snr.R
# ==============================================================================
#
# Analytic data generating process mirroring the hapgen_like.py LD block
# correlation structure. Sourced by demo files — do not run standalone.
#
# Block structure (1-based R column indices, p = 150 SNPs):
#
#   Block 1:  SNPs   1– 45  base_corr = 0.40  (moderate LD)
#   Block 2:  SNPs  46– 80  base_corr = 0.80  (strong LD)
#   Block 3:  SNPs  81– 90  base_corr = 0.95  (very strong LD)
#   Block 4:  SNPs  91–110  base_corr = 0.00  (independent)
#   Block 5:  SNPs 111–140  base_corr = 0.70  (strong LD)
#   Block 6:  SNPs 141–150  base_corr = 0.50  (moderate-strong LD)
#
# Within-block correlation model (analytic, no random signs / noise):
#
#   Sigma[i, j] = base_corr * exp(-|rank_i - rank_j| * decay_rate)
#   Sigma[i, i] = 1   (unit diagonal)
#
# where rank_i and rank_j are 1-based positions within the block.
# This is equivalent to an AR(1)-type exponential decay with
# phi = exp(-decay_rate), matching the formula in hapgen_like.py
# (np.exp(-dist * 0.1)) but without the random per-pair sign or additive
# noise.  Long-range LD from hapgen_like.py (random 50/50 sign, blocks
# (10-19, 120-129) and (50-59, 90-99)) has E[cor] ≈ 0 and is omitted in
# this analytic version; between-block entries are zero.
#
# Active variables: the first s_per_block column indices of each block.
#
# Public interface:
#   make_hapgen_blocks()    — named list of 6 block definitions
#   make_hapgen_sigma()     — 150 x 150 block correlation matrix
#   make_hapgen_support()   — active predictor index vector (1-based)
#   dgp_hapgen_snr()        — full DGP: returns list(X, y, beta, ...)
# ==============================================================================

# Fixed number of SNPs (matches hapgen_like.py)
.HAPGEN_P <- 150L


# ==============================================================================
# make_hapgen_blocks
# ==============================================================================

#' Return the 6 LD block definitions used in the hapgen-like DGP.
#'
#' Each element is a named list with fields:
#'   start     — first 1-based column index of the block
#'   end       — last  1-based column index of the block
#'   base_corr — within-block base correlation (0 = independent block)
#'
#' Block proportions from hapgen_like.py (template: p0 = 150 SNPs):
#'   Block 1: 45/150 = 30%   base_corr = 0.40  (moderate LD)
#'   Block 2: 35/150 ≈ 23%   base_corr = 0.80  (strong LD)
#'   Block 3: 10/150 ≈  7%   base_corr = 0.95  (very strong LD)
#'   Block 4: 20/150 ≈ 13%   base_corr = 0.00  (independent)
#'   Block 5: 30/150 = 20%   base_corr = 0.70  (strong LD)
#'   Block 6: 10/150 ≈  7%   base_corr = 0.50  (moderate-strong LD)
#' When p != 150, block sizes are scaled via round(raw / 150 * p) with
#' any rounding residual absorbed by the last block.
#'
#' @param p  Total number of SNPs. Default 150.
#' @return List of 6 named lists.
make_hapgen_blocks <- function(p = 150L) {
  raw_sizes  <- c(45L, 35L, 10L, 20L, 30L, 10L)   # proportions at p0 = 150
  base_corrs <- c(0.40, 0.80, 0.95, 0.00, 0.70, 0.50)

  if (p == 150L) {
    sizes <- raw_sizes
  } else {
    sizes     <- as.integer(round(raw_sizes / 150.0 * p))
    sizes[6L] <- sizes[6L] + (p - sum(sizes))   # absorb rounding residual
  }

  ends   <- cumsum(sizes)
  starts <- c(1L, ends[-6L] + 1L)

  lapply(seq_len(6L), function(i) {
    list(start = starts[i], end = ends[i], base_corr = base_corrs[i])
  })
}


# ==============================================================================
# make_hapgen_sigma
# ==============================================================================

#' Build the 150 x 150 analytic hapgen-like LD block correlation matrix.
#'
#' Within each non-independent block m the (i, j) entry is:
#'
#'   Sigma[i, j] = base_corr_m * exp(-|rank(i,m) - rank(j,m)| * decay_rate)
#'
#' where rank(i, m) is the 1-based position of SNP i within block m.
#' The diagonal is 1.  All between-block entries are 0 (block diagonal).
#'
#' Positive definiteness is guaranteed analytically: within each block the
#' matrix equals (1 - c) * I + c * Phi where c = base_corr and Phi is the
#' AR(1) Toeplitz matrix with phi = exp(-decay_rate) in (0, 1), which is
#' positive definite.  A small jitter is added to the diagonal for numerical
#' robustness.
#'
#' @param decay_rate  Exponential decay rate per unit rank distance. Default
#'   0.1, matching hapgen_like.py (np.exp(-dist * 0.1)).
#' @param jitter      Scalar added to every diagonal entry after block
#'   construction. Default 1e-6.
#'
#' @param p  Total number of SNPs. Default 150.
#' @return p x p numeric matrix (symmetric, positive definite).
#'
#' @examples
#' Sigma <- make_hapgen_sigma()            # 150 x 150 (default)
#' Sigma <- make_hapgen_sigma(p = 1000L)  # 1000 x 1000
#' min(eigen(Sigma, only.values = TRUE)$values)  # > 0
make_hapgen_sigma <- function(p = 150L, decay_rate = 0.1, jitter = 1e-6) {

  Sigma  <- diag(p)
  blocks <- make_hapgen_blocks(p)

  for (bl in blocks) {
    if (bl$base_corr == 0) next

    idx <- seq.int(bl$start, bl$end)
    b_m <- length(idx)
    r   <- seq_len(b_m)

    within_block <- outer(r, r, function(i, j) {
      bl$base_corr * exp(-abs(i - j) * decay_rate)
    })
    diag(within_block) <- 1.0

    Sigma[idx, idx] <- within_block
  }

  # Numerical PD guarantee
  diag(Sigma) <- diag(Sigma) + jitter
  Sigma
}


# ==============================================================================
# make_hapgen_support
# ==============================================================================

#' Build the active predictor support for the hapgen-like DGP.
#'
#' Selects the first s_per_block column indices (1-based) from each of the
#' 6 LD blocks, yielding s_total = 6 * s_per_block active predictors.
#' This mirrors the "first s variables per group" convention in dgp_gvs.R.
#'
#' @param s_per_block  Number of active predictors per LD block. Must satisfy
#'   1 <= s_per_block <= min block size = 10 (block 3 and block 6).
#'   Default 1L.
#'
#' @return Sorted integer vector of length 6 * s_per_block (1-based).
#'
#' @examples
#' make_hapgen_support(1L)
#' # -> c(1, 46, 81, 91, 111, 141)  (first predictor of each block)
#'
#' make_hapgen_support(2L)
#' # -> c(1, 2, 46, 47, 81, 82, 91, 92, 111, 112, 141, 142)
make_hapgen_support <- function(s_per_block = 1L, p = 150L) {
  blocks <- make_hapgen_blocks(p)
  stopifnot(s_per_block >= 1L)

  min_block_size <- min(vapply(blocks, function(bl) bl$end - bl$start + 1L, integer(1L)))
  if (s_per_block > min_block_size)
    stop(sprintf("'s_per_block' (%d) exceeds minimum block size (%d).",
                 s_per_block, min_block_size))

  support <- unlist(lapply(blocks, function(bl) {
    seq.int(bl$start, bl$start + s_per_block - 1L)
  }))
  sort(as.integer(support))
}


# ==============================================================================
# dgp_hapgen_snr
# ==============================================================================

#' DGP for the hapgen-like LD block scenario with SNR-calibrated noise.
#'
#' Predictor rows are drawn i.i.d. from N(0, Sigma) where Sigma is the
#' 150 x 150 analytic LD block correlation matrix built by make_hapgen_sigma().
#' The matrix is drawn via the lower Cholesky factor X = W L' and
#' column-centred to remove finite-sample mean shifts.
#'
#' The response follows:
#'
#'   y = X %*% beta + eps,   eps ~ N(0, noise_sigma^2 * I_n)
#'
#' where beta[support] = amplitude and noise_sigma is chosen so that
#'
#'   Var(X %*% beta) / Var(eps) = snr
#'
#' using the n-1 denominator for both variances, matching dgp_ar1_snr().
#'
#' @param n          Number of observations.
#' @param p          Total number of SNPs (predictors). Default 150 (matches
#'                   hapgen_like.py). The 6-block LD structure is scaled
#'                   proportionally to p via make_hapgen_blocks(p).
#' @param support    1-based integer vector of active predictor indices.
#'                   Typically built via make_support_capped_spread(s, p, max_gap).
#' @param amplitude  Signal coefficient for every active predictor. Default 3.
#' @param snr        Target signal-to-noise ratio Var(signal) / Var(noise).
#'                   Default 2.
#' @param decay_rate Exponential decay rate for within-block correlation.
#'                   Default 0.1 (matches hapgen_like.py).
#' @param jitter     Diagonal regularisation added to Sigma. Default 1e-6.
#' @param seed       Integer RNG seed. Default NULL.
#'
#' @return Named list:
#'   \describe{
#'     \item{X}{Predictor matrix (n x p), column-centred.}
#'     \item{y}{Response vector (length n).}
#'     \item{beta}{Coefficient vector (length p; non-zero at support.)}
#'     \item{true_support}{1-based active predictor indices (integer vector).}
#'     \item{n, p, s, snr}{Dimension and parameter scalars.}
#'     \item{Sigma}{p x p population covariance used to draw X.}
#'   }
#'
#' @examples
#' dat <- dgp_hapgen_snr(
#'   n         = 300,
#'   p         = 1000L,
#'   support   = make_support_capped_spread(10L, 1000L, 5L),
#'   amplitude = 2,
#'   snr       = 2.0,
#'   seed      = 2026
#' )
#' dim(dat$X)          # 300 x 1000
#' dat$true_support    # 1 6 11 16 21 26 31 36 41 46
dgp_hapgen_snr <- function(n,
                           p          = 150L,
                           support,
                           amplitude  = 3.0,
                           snr        = 2.0,
                           decay_rate = 0.1,
                           jitter     = 1e-6,
                           seed       = NULL) {

  stopifnot(is.numeric(support) || is.integer(support))
  stopifnot(all(support >= 1L), all(support <= p))
  stopifnot(snr > 0, amplitude > 0)

  if (!is.null(seed)) set.seed(seed)

  # ── Covariance matrix and Cholesky draw ───────────────────────────────────────
  Sigma <- make_hapgen_sigma(p = p, decay_rate = decay_rate, jitter = jitter)
  L     <- t(chol(Sigma))                              # lower Cholesky factor
  W     <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X     <- W %*% t(L)
  X     <- X - matrix(colMeans(X), nrow = n, ncol = p, byrow = TRUE)  # centre

  # ── Coefficient vector ────────────────────────────────────────────────────────
  beta          <- numeric(p)
  beta[support] <- amplitude

  # ── SNR-calibrated response ───────────────────────────────────────────────────
  # noise_sigma = sqrt(var(signal) / snr)   (n-1 denominator, Bessel corrected)
  signal      <- drop(X %*% beta)
  var_sig     <- stats::var(signal)
  noise_sigma <- sqrt(var_sig / snr)
  y           <- signal + stats::rnorm(n, sd = noise_sigma)

  list(
    X            = X,
    y            = y,
    beta         = beta,
    true_support = as.integer(support),
    n            = n,
    p            = p,
    s            = length(support),
    snr          = snr,
    Sigma        = Sigma
  )
}
