# ==============================================================================
# dgp_arp_snr.R
#
# Data generating processes for the AR(p) scenario, using SNR-based noise
# calibration consistent with dgp_ar1_snr.R and da_trex_sim_common.hpp.
#
# The AR(1) special case reproduces dgp_ar1_snr() exactly (same RNG sequence,
# same column-by-column recursion, same SNR calibration via sample var).
#
# AR(p) model (p >= 1):
#
#   X[i, j] = phi[1]*X[i,j-1] + phi[2]*X[i,j-2] + ... + phi[p]*X[i,j-p]
#             + sigma_eps * eta[i, j]
#
#   with eta[i, j] i.i.d. N(0, 1) and sigma_eps chosen so that Var(X[i,j]) = 1
#   (unit-variance columns) via the Yule-Walker equation:
#
#   sigma_eps^2 = 1 - phi %*% r[1..p]
#
#   where r[d] = Corr(X[i,j], X[i,j-d]) is the theoretical ACF, computed
#   by solving the Yule-Walker linear system (see .ar_acf_yw()).
#
# Coefficient families provided:
#
#   make_ar_phi_geometric(p, rho, decay)
#     phi[l] = rho * decay^(l-1),  l = 1..p
#     All share phi[1] = rho; higher lags decay geometrically.
#     Stationarity is checked; an error is raised if violated.
#
#   make_ar_phi_pacf(p, rho)
#     Levinson-Durbin from constant PACF = rho at all lags.
#     Guarantees stationarity for any rho in (-1, 1) and any p.
#
# Public functions:
#   make_ar_phi_geometric(p, rho, decay)
#   make_ar_phi_pacf(p, rho)
#   dgp_arp_snr(n, p, support, amplitude, phi, snr, seed)
#
# Relationship to dgp_ar1_snr():
#   dgp_arp_snr(..., phi = c(rho)) reproduces dgp_ar1_snr(..., rho = rho)
#   column-for-column and RNG-call-for-RNG-call.
#
# ==============================================================================


# ==============================================================================
# Internal helpers
# ==============================================================================

#' Theoretical ACF r[0..max_lag] for a stationary AR(phi) process.
#'
#' Solves the Yule-Walker linear system for r[1..p] exactly, then recurses
#' for higher lags.  Unit variance is assumed (r[0] = 1).
#'
#' @param phi Numeric vector of AR coefficients phi[1..p].
#' @param max_lag Maximum lag to compute.
#' @return Numeric vector of length max_lag + 1 with r[0], r[1], ..., r[max_lag].
.ar_acf_yw <- function(phi, max_lag) {
  p <- length(phi)

  # Build the (p x p) Yule-Walker system for r[1..p]:
  #   r[k] = sum_{l=1}^p  phi[l] * r[|k - l|],   k = 1..p,  r[0] = 1
  # Rearranged into A * [r1..rp]^T = b:
  A <- matrix(0.0, nrow = p, ncol = p)
  b <- numeric(p)
  for (k in seq_len(p)) {
    b[k] <- phi[k]                        # phi[k] * r[0] term
    for (l in seq_len(p)) {
      lag <- abs(k - l)
      if (lag == 0L) {
        A[k, l] <- A[k, l] + 1.0         # r[k] on LHS
      } else {
        A[k, lag] <- A[k, lag] - phi[l]  # move phi[l]*r[lag] to LHS
      }
    }
  }
  r_init <- solve(A, b)

  r <- numeric(max_lag + 1L)
  r[1L] <- 1.0                            # r[0]
  for (k in seq_len(min(p, max_lag)))
    r[k + 1L] <- r_init[k]               # r[1..p]
  for (k in seq(p + 1L, max_lag)) {
    r[k + 1L] <- sum(phi * r[k - seq_len(p) + 1L])
  }
  r
}


#' Check stationarity of AR(phi): all roots of the AR polynomial outside unit circle.
#'
#' @param phi Numeric vector of AR coefficients.
#' @return Maximum modulus of the companion-matrix eigenvalues.
.ar_max_root <- function(phi) {
  p <- length(phi)
  if (p == 1L) return(abs(phi[1L]))
  C <- matrix(0.0, nrow = p, ncol = p)
  C[1L, ] <- phi
  C[2L:p, 1L:(p - 1L)] <- diag(p - 1L)
  max(Mod(eigen(C, only.values = TRUE)$values))
}


# ==============================================================================
# Coefficient constructors
# ==============================================================================

#' Geometric AR(p) coefficients: phi[l] = rho * decay^(l-1).
#'
#' All process orders share phi[1] = rho, with higher lags decaying at rate
#' `decay`.  The function checks stationarity and stops if violated.
#'
#' Suitable decay values for rho = 0.7:
#'   p = 2: decay <= 0.40
#'   p = 3: decay <= 0.30
#'   p = 5: decay <= 0.30
#' (tighter decay required for higher orders)
#'
#' @param p    AR order (>= 1).
#' @param rho  Lag-1 coefficient phi[1] = rho.  In (0, 1).
#' @param decay Geometric decay rate for phi[l] = rho * decay^(l-1). Default 0.3.
#'
#' @return Named numeric vector phi[1..p].
#'
#' @examples
#' make_ar_phi_geometric(1, rho = 0.7)          # -> c(phi1 = 0.7)
#' make_ar_phi_geometric(2, rho = 0.7, decay = 0.3)
#' make_ar_phi_geometric(5, rho = 0.7, decay = 0.3)
make_ar_phi_geometric <- function(p, rho, decay = 0.3) {
  stopifnot(p >= 1L)
  stopifnot(rho > 0, rho < 1)
  stopifnot(decay > 0, decay < 1)
  phi <- rho * decay^(seq_len(p) - 1L)
  names(phi) <- paste0("phi", seq_len(p))
  mr <- .ar_max_root(phi)
  if (mr >= 1.0)
    stop(sprintf(
      "AR(%d) with rho=%.2f, decay=%.2f is non-stationary (max|root|=%.4f). ",
      p, rho, decay,  mr,
      "Reduce `decay` or `rho`."
    ))
  phi
}


#' Levinson-Durbin AR(p) coefficients from constant PACF = rho.
#'
#' Sets the partial autocorrelation at every lag equal to rho, then converts
#' to AR coefficients via the Levinson-Durbin recursion.  Stationarity is
#' guaranteed for any rho in (-1, 1).
#'
#' Note: for p = 1 this reduces to phi = c(rho), identical to AR(1).
#' For p > 1 the lag-1 ACF exceeds rho (it is amplified by the higher-order
#' terms), so this family does NOT share phi[1] = rho across orders.
#'
#' @param p   AR order (>= 1).
#' @param rho Constant PACF value at all lags.  In (-1, 1).
#'
#' @return Named numeric vector phi[1..p].
#'
#' @examples
#' make_ar_phi_pacf(1, rho = 0.7)   # -> c(phi1 = 0.7)
#' make_ar_phi_pacf(2, rho = 0.7)
#' make_ar_phi_pacf(5, rho = 0.5)
make_ar_phi_pacf <- function(p, rho) {
  stopifnot(p >= 1L)
  stopifnot(abs(rho) < 1)
  phi_mat <- matrix(0.0, nrow = p, ncol = p)
  phi_mat[1L, 1L] <- rho
  for (k in seq(2L, p)) {
    pkk <- rho
    phi_mat[k, k] <- pkk
    for (j in seq_len(k - 1L))
      phi_mat[k, j] <- phi_mat[k - 1L, j] - pkk * phi_mat[k - 1L, k - j]
  }
  phi <- phi_mat[p, ]
  names(phi) <- paste0("phi", seq_len(p))
  phi
}


# ==============================================================================
# dgp_arp_snr
# ==============================================================================

#' AR(p) DGP with SNR-calibrated noise.
#'
#' Generates a predictor matrix X (n x p) whose columns follow a stationary
#' AR(p) process with coefficient vector `phi`, then builds an SNR-calibrated
#' linear response y.
#'
#' Column generation (mirrors dgp_ar1_snr for p = 1):
#'
#'   X[i, 1..q]  drawn i.i.d. N(0, 1)          (initialisation, q = length(phi))
#'   X[i, j]   = sum_{l=1}^q phi[l]*X[i,j-l]
#'               + sigma_eps * eta[i, j],        j = q+1..p
#'
#' with sigma_eps = sqrt(1 - phi %*% r[1..q]) derived from the theoretical ACF
#' via Yule-Walker, so that Var(X[i,j]) = 1 in the stationary regime.
#'
#' The first q columns are a short boundary artefact (drawn from the marginal
#' N(0,1) rather than the joint stationary distribution).  For typical values
#' p >> q this is negligible.  Use `burn_in > 0` to discard the first
#' `burn_in` columns and replace them with stationary draws if needed.
#'
#' Response:
#'
#'   y = X %*% beta + eps,   eps ~ N(0, sigma_noise^2 * I_n)
#'   sigma_noise = sqrt( sample_var(X %*% beta) / snr )
#'
#' with sample_var() using the n-1 denominator, consistent with
#' dgp_ar1_snr() and the C++ make_snr_response().
#'
#' @param n         Number of observations.
#' @param p         Number of predictors.
#' @param support   1-based integer vector of active predictor indices.
#' @param amplitude Signal coefficient for each active predictor. Default 3.
#' @param phi       AR coefficient vector phi[1..q] (length q = AR order).
#'                  Use make_ar_phi_geometric() or make_ar_phi_pacf() to build.
#'                  Default c(0.7) reproduces dgp_ar1_snr() with rho = 0.7.
#' @param snr       Target SNR = Var(signal) / Var(noise). Default 2.
#' @param seed      Integer RNG seed (governs both X and noise). Default NULL.
#' @param burn_in   Number of leading columns to discard for stationarity.
#'                  Default 0 (consistent with dgp_ar1_snr convention).
#'
#' @return A list with components:
#'   X            Predictor matrix (n x p).
#'   y            Response vector (length n).
#'   beta         Coefficient vector (length p).
#'   true_support 1-based integer index vector of active predictors.
#'   n, p, s      Dimensions.
#'   phi          AR coefficient vector used.
#'   ar_order     AR order q = length(phi).
#'   snr          Target SNR.
#'   sigma_eps    Innovation standard deviation.
#'   acf_theory   Theoretical ACF r[0..q] (length q+1).
#'
#' @examples
#' # AR(1) — identical to dgp_ar1_snr(rho = 0.7):
#' support <- make_support_capped_spread(s = 10L, p = 1000L, max_gap = 20L)
#' dat1 <- dgp_arp_snr(n = 300, p = 1000, support = support,
#'                     phi = make_ar_phi_geometric(1, rho = 0.7),
#'                     snr = 2.0, seed = 2026)
#'
#' # AR(2):
#' dat2 <- dgp_arp_snr(n = 300, p = 1000, support = support,
#'                     phi = make_ar_phi_geometric(2, rho = 0.7, decay = 0.3),
#'                     snr = 2.0, seed = 2026)
#'
#' # AR(5) via PACF parametrisation:
#' dat5 <- dgp_arp_snr(n = 300, p = 1000, support = support,
#'                     phi = make_ar_phi_pacf(5, rho = 0.5),
#'                     snr = 2.0, seed = 2026)
dgp_arp_snr <- function(n,
                        p,
                        support,
                        amplitude = 3,
                        phi       = c(0.7),
                        snr       = 2.0,
                        seed      = NULL,
                        burn_in   = 0L) {

  stopifnot(is.integer(support) || is.numeric(support))
  stopifnot(all(support >= 1L), all(support <= p))
  stopifnot(snr > 0)
  stopifnot(burn_in >= 0L)

  q <- length(phi)

  # Stationarity check
  mr <- .ar_max_root(phi)
  if (mr >= 1.0)
    stop(sprintf("AR(%d) is non-stationary (max|root| = %.4f). ", q, mr,
                 "Use make_ar_phi_geometric() or make_ar_phi_pacf() to ensure stationarity."))

  # Theoretical ACF and innovation std
  acf_theory <- .ar_acf_yw(phi, max_lag = q)
  sigma_eps   <- sqrt(max(1.0 - sum(phi * acf_theory[2L:(q + 1L)]), 0.0))

  if (!is.null(seed)) set.seed(seed)

  # Generate X: AR(q) recursion column-by-column
  # RNG call order: col 1, col 2, ..., col (p + burn_in)
  # For q = 1, phi = c(rho): reproduces dgp_ar1_snr() exactly.
  # ------------------------------------------------------------------
  p_total <- p + burn_in
  X_full  <- matrix(0.0, nrow = n, ncol = p_total)

  # Initialise first q columns from N(0, 1)
  for (j in seq_len(q))
    X_full[, j] <- stats::rnorm(n)

  # Stationary recursion for columns q+1 .. p_total
  for (j in seq(q + 1L, p_total)) {
    col_j <- sigma_eps * stats::rnorm(n)
    for (l in seq_len(q))
      col_j <- col_j + phi[l] * X_full[, j - l]
    X_full[, j] <- col_j
  }

  # Discard burn-in columns
  X <- X_full[, (burn_in + 1L):p_total, drop = FALSE]

  # Coefficient vector
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
    s            = length(support),
    phi          = phi,
    ar_order     = q,
    snr          = snr,
    sigma_eps    = sigma_eps,
    acf_theory   = acf_theory
  )
}
