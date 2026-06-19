# ==============================================================================
# dgp_hapgen_snr.R
# ==============================================================================
# Data-generating process for the quasi-hapgen LD-block scenario.
# Sourced by demo_trex_gvs_12.R.
#
# Contents:
#   dgp_hapgen_snr(n, snr, rho_scale, seed)
#       SNR-calibrated DGP — 7 irregular LD blocks (heterogeneous AR(1)-like
#       within-block correlation) + two weak long-range cross-block LD pairs.
#       Fixed p = 500. Active set: blocks 1, 3, 4 (s = 130).
#
# Block layout (p = 500, no shuffling — fixed genomic-style order):
#
#   Block 1: cols   1- 30  size  30  rho_base=0.85  ACTIVE
#   gap:     cols  31- 35  size   5  white noise
#   Block 2: cols  36-115  size  80  rho_base=0.45  trap
#   gap:     cols 116-120  size   5  white noise
#   Block 3: cols 121-180  size  60  rho_base=0.70  ACTIVE
#   gap:     cols 181-185  size   5  white noise
#   Block 4: cols 186-225  size  40  rho_base=0.90  ACTIVE
#   gap:     cols 226-230  size   5  white noise
#   Block 5: cols 231-330  size 100  rho_base=0.55  trap
#   gap:     cols 331-335  size   5  white noise
#   Block 6: cols 336-385  size  50  rho_base=0.80  trap
#   gap:     cols 386-390  size   5  white noise
#   Block 7: cols 391-500  size 110  rho_base=0.30  trap (background)
#
# Long-range LD pairs (scaled by rho_scale):
#   cols  10-15  (Block 1)  <->  cols 340-345  (Block 6): strength 0.30
#   cols 130-135 (Block 3)  <->  cols 240-245  (Block 5): strength 0.25
#
# Sigma construction:
#   1. I_500 base
#   2. Each block: Sigma[b, b] = rho_eff_b^|i-j|, rho_eff_b = min(rho_base * rho_scale, 0.999)
#   3. Inject long-range LD (also scaled by rho_scale, clipped to 0.999)
#   4. Symmetrize; repair PSD via Matrix::nearPD(corr = TRUE)
# ==============================================================================


# ==============================================================================
# dgp_hapgen_snr
# ==============================================================================
#' SNR-calibrated quasi-hapgen LD-block DGP.
#'
#' Fixed p = 500. 7 irregular LD blocks (heterogeneous AR(1)-like correlation) +
#' two weak long-range cross-block LD pairs. Active set: blocks 1, 3, 4 (s = 130).
#'
#' @param n         Number of observations.
#' @param snr       Signal-to-noise ratio (var(signal) / var(noise)).
#' @param rho_scale Global correlation scale in [0, 1]. Multiplies all within-block
#'                  base correlations and long-range LD strengths. Default 1.0.
#' @param seed      RNG seed. Default NULL.
#' @return list with fields:
#'   X, y, beta, n, p, s, snr, sigma_y, rho_scale,
#'   block_indices, block_sizes, block_rhos_base, block_rhos_eff
dgp_hapgen_snr <- function(n, snr, rho_scale = 1.0, seed = NULL) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for dgp_hapgen_snr().")
  }

  if (!is.null(seed)) set.seed(seed)

  # ============================================================
  # Fixed layout (no shuffling)
  # ============================================================
  p                <- 500L
  block_starts     <- c(  1L,  36L, 121L, 186L, 231L, 336L, 391L)
  block_ends       <- c( 30L, 115L, 180L, 225L, 330L, 385L, 500L)
  block_sizes      <- block_ends - block_starts + 1L
  block_rhos_base  <- c(0.85, 0.45, 0.70, 0.90, 0.55, 0.80, 0.30)
  is_active        <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  n_blocks         <- length(block_starts)

  block_rhos_eff <- pmin(block_rhos_base * rho_scale, 0.999)

  block_indices <- vector("list", n_blocks)
  for (b in seq_len(n_blocks)) {
    block_indices[[b]] <- block_starts[b]:block_ends[b]
  }

  # ============================================================
  # Build Sigma
  # ============================================================
  Sigma <- diag(p)

  for (b in seq_len(n_blocks)) {
    idx    <- block_indices[[b]]
    rho_b  <- block_rhos_eff[b]
    sz     <- block_sizes[b]
    r_idx  <- seq_len(sz)
    Sigma_b <- rho_b^abs(outer(r_idx, r_idx, "-"))
    Sigma[idx, idx] <- Sigma_b
  }

  # Long-range LD injection
  lr_pairs <- list(
    list(rows = 10L:15L, cols = 340L:345L, strength = 0.30),
    list(rows = 130L:135L, cols = 240L:245L, strength = 0.25)
  )
  for (lr in lr_pairs) {
    val <- min(lr$strength * rho_scale, 0.999)
    Sigma[lr$rows, lr$cols] <- val
    Sigma[lr$cols, lr$rows] <- val
  }

  # Symmetrize and repair PSD
  Sigma  <- (Sigma + t(Sigma)) / 2
  Sigma_psd <- as.matrix(Matrix::nearPD(Sigma, corr = TRUE)$mat)

  # ============================================================
  # Simulate X
  # ============================================================
  chol_S <- chol(Sigma_psd)           # upper-triangular R s.t. R'R = Sigma_psd
  Z      <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X      <- Z %*% chol_S              # rows ~ MVN(0, Sigma_psd)

  # ============================================================
  # Beta and response
  # ============================================================
  beta <- rep(0, p)
  for (b in seq_len(n_blocks)) {
    if (is_active[b]) beta[block_indices[[b]]] <- 3
  }

  signal  <- as.numeric(X %*% beta)
  sigma_y <- sqrt(var(signal) / snr)
  y       <- signal + rnorm(n, sd = sigma_y)

  list(
    X               = X,
    y               = y,
    beta            = beta,
    n               = n,
    p               = p,
    s               = 130L,
    snr             = snr,
    sigma_y         = sigma_y,
    rho_scale       = rho_scale,
    block_indices   = block_indices,
    block_sizes     = block_sizes,
    block_rhos_base = block_rhos_base,
    block_rhos_eff  = block_rhos_eff
  )
}
