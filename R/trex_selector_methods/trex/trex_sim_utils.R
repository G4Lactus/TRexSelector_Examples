# ==============================================================================
# trex_sim_utils.R
# ==============================================================================
#
# Shared simulation utilities for T-Rex demos.
#
# Mirrors C++ trex_sim_utils.hpp.
#
# Contents:
# -----------------------------------------------------------------------------
#   Section 1 — DGP (in-memory)
#     dgp_gauss_snr()           — i.i.d. Gaussian DGP with SNR-calibrated noise.
#
#   Section 2 — DGP (out-of-core / mmap)
#     dgp_chunked_gauss_snr()   — Chunked out-of-core DGP for massive datasets.
#
#   Section 3 — Single-run print helpers
#     .print_single_run()       — Formatted single-run result block.
#
#   Section 4 — MC infrastructure
#     SOLVERS_DEFAULT            — List of 14 solver descriptors (mirrors
#                                  make_default_solvers_to_test() in C++).
#     make_mmap_trex_ctrl()      — Build trex_control for memory-mapped scenarios.
#     execute_with_temp_mmap()   — RAII-like lifecycle manager for temp mmap files.
#     .print_demo_header()       — Standardised demo section headers.
#     .make_trex_ctrl()          — Build trex_control() from a solver descriptor.
#     .run_mc_trex()             — Parallel MC runner (TRexSelector OOP API).
#     .print_solver_table()      — Aligned solver x SNR table to console.
#     .save_mc_csv()             — Tidy long-format CSV to disk.
#     .save_and_print_mc()       — Console + .txt + .csv combined output.
#
# ==============================================================================


# ==============================================================================
# Section 1 — DGP (in-memory)
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
  signal      <- drop(X[, support] %*% beta[support])
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


# ==============================================================================
# Section 2 — DGP (out-of-core / mmap)
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
  X_mmap <- TRexSelectorNeo::mmap_matrix(x_path, n, p)

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


# ==============================================================================
# Section 3 — Single-run print helpers
# ==============================================================================

#' Print a formatted single-run result block.
.print_single_run <- function(scenario_name, dat, selector, tFDR) {

  true_support <- dat$true_support
  s            <- length(true_support)
  p            <- ncol(dat$X)
  selected     <- selector$selected_indices
  n_sel        <- length(selected)
  n_tp         <- length(intersect(selected, true_support))
  n_fp         <- n_sel - n_tp
  tpp          <- TRexSelectorNeo::compute_tpp(selected, true_support)
  fdp          <- TRexSelectorNeo::compute_fdp(selected, true_support)

  cat(strrep("=", 70L), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70L), "\n")
  cat(sprintf("  Data:  n = %d,  p = %d,  s = %d,  tFDR = %.2f\n",
              nrow(dat$X), p, s, tFDR))
  cat(sprintf("         True support (1-based):   {%s}\n",
              paste(true_support, collapse = ", ")))
  cat(sprintf("         Active coefficients:       {%s}\n",
              paste(round(dat$beta[true_support], 3L), collapse = ", ")))
  cat(strrep("-", 70L), "\n")
  cat(sprintf("  Calibration:  T_stop = %d,  L = %d\n",
              selector$T_stop, selector$L))
  cat(sprintf("  Selection:    %d selected  |  TP = %d  FP = %d\n",
              n_sel, n_tp, n_fp))
  cat(sprintf("  Rates:        TPP = %.4f  |  FDP = %.4f  (target tFDR <= %.2f)\n",
              tpp, fdp, tFDR))
  cat(strrep("-", 70L), "\n")
  cat("\nAdjusted Relative Occurrences (phi_prime):\n")
  cat(format(round(selector$phi_prime, 6L)), fill = TRUE)
  cat("\nPhi Matrix (phi_mat):\n")
  print(round(selector$phi_mat, 4L))
  cat("\nEstimated FDP Matrix (fdp_hat_mat):\n")
  print(round(selector$fdp_hat_mat, 4L))
  cat("\nR Matrix (r_mat):\n")
  print(round(selector$r_mat, 4L))
  cat("\nVoting Grid (voting_grid):\n")
  cat(format(round(selector$voting_grid, 4L)), fill = TRUE)
  cat(strrep("=", 70L), "\n\n")
  invisible(NULL)

}


# ==============================================================================
# Section 4 — MC infrastructure
# ==============================================================================
#
# Shared infrastructure for classic T-Rex MC simulations.
# Sourced by demo files — do not run standalone.
#
# Mirrors C++ trex_sim_utils.hpp.


# ==============================================================================
# Default solver list
# (mirrors make_default_solvers_to_test() in C++)
# ==============================================================================

SOLVERS_DEFAULT <- list(
  list(solver = "TLARS",      name = "TLARS",        lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TLASSO",     name = "TLASSO",       lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TENET",      name = "TENET",        lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TSTAGEWISE", name = "TSTAGEWISE",   lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TSTEPWISE",  name = "TSTEPWISE",    lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TOMP",       name = "TOMP",         lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TGP",        name = "TGP",          lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TACGP",      name = "TACGP",        lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TMP",        name = "TMP",          lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TAFS",       name = "TAFS_rho_0.3", lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TAFS",       name = "TAFS_rho_1.0", lambda2 = 0.1, rho_afs = 1.0,
       ncgmp_variant = 0L),
  list(solver = "TNCGMP",     name = "TNCGMP_v1",    lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 1L),
  list(solver = "TNCGMP",     name = "TNCGMP_v0",    lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L),
  list(solver = "TOOLS",      name = "TOOLS",        lambda2 = 0.1, rho_afs = 0.3,
       ncgmp_variant = 0L)
)


# ==============================================================================
# Shared Helpers for Memory Mapping and Console Output
# ==============================================================================

#' Make a TRex control list for memory-mapping scenarios.
#'
#' Matches C++ make_mmap_trex_ctrl().
#'
#' @param ... Optional overrides for trex_control parameters.
#' @return A trex_control object.
make_mmap_trex_ctrl <- function(...) {

  base_ctrl <- list(
    K                        = 20L,
    max_dummy_multiplier     = 10L,
    use_max_T_stop           = TRUE,
    lloop_strategy           = "HCONCAT",
    tloop_stagnation_stop    = TRUE,
    use_memory_mapping       = TRUE
  )

  overrides <- list(...)
  for (nm in names(overrides)) {
    base_ctrl[[nm]] <- overrides[[nm]]
  }
  do.call(TRexSelectorNeo::trex_control, base_ctrl)
}


#' Execute a block of code with a memory-mapped temporary file
#'
#' Handles the RAII-like lifecycle of writing an in-memory matrix to a
#' temporary memory-mapped file, yielding it to the provided function,
#' and ensuring the file is unlinked safely.
#'
#' @param X_matrix The in-memory matrix to convert.
#' @param expr_func A function taking the memory-mapped matrix as its argument.
#' @return The result of `expr_func`.
execute_with_temp_mmap <- function(X_matrix, expr_func) {

  x_path <- tempfile(fileext = ".dat")
  on.exit(unlink(x_path), add = TRUE)

  X_mmap <- TRexSelectorNeo::convert_to_memory_mapped(X_matrix, x_path)
  expr_func(X_mmap)

}


#' Print standardized demo section headers
#'
#' @param part_name Section title.
#' @param description Short description.
#' @param high_dim Boolean flag for high vs low dimensional setup.
#' @param n Number of samples.
#' @param p Number of features.
#' @param num_MC Number of Monte Carlo trials (optional).
.print_demo_header <- function(part_name, description, high_dim, n, p, num_MC = NULL) {

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"

  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat(part_name, "\n")
  cat(sprintf("  %s\n", description))

  if (is.null(num_MC)) {
    cat(sprintf("  %s  (n = %d, p = %d)\n\n", dim_label, n, p))
  } else {
    cat(sprintf("  %s  (n = %d, p = %d,  num_MC = %d)\n\n",
                dim_label, n, p, num_MC))
  }
}


# ==============================================================================
# .make_trex_ctrl
# ==============================================================================

#' Build a trex_control() object from a solver descriptor and a base-params list.
#'
#' @param solver    One element of SOLVERS_DEFAULT (a named list with solver,
#'                  lambda2, rho_afs, ncgmp_variant).
#'
#' @param base      Named list of shared control parameters forwarded verbatim
#'                  to trex_control(). Must NOT include solver, lambda2,
#'                  rho_afs or ncgmp_variant (those come from `solver`).
#'
#' @return A trex_control() object.
.make_trex_ctrl <- function(solver, base) {
  do.call(
    trex_control,
    c(list(solver        = solver$solver,
           lambda2       = solver$lambda2,
           rho_afs       = solver$rho_afs,
           ncgmp_variant = solver$ncgmp_variant),
      base)
  )
}


# ==============================================================================
# .run_mc_trex
# ==============================================================================

#' Parallel MC runner for classical T-Rex (TRexSelector OOP API).
#'
#' Runs `num_MC` independent trials. Each trial calls `make_data_fn(mc,
#' trial_seed)` to obtain the data, builds a TRexSelector with `ctrl`, and
#' records FDP / TPP (and optionally L and T_stop).
#'
#' Parallelism uses `parallel::mclapply` on Unix (fork-based; each child has
#' its own R process so set.seed() inside make_data_fn is thread-safe) and
#' `parallel::parLapply` on Windows.
#'
#' @param make_data_fn  Factory: function(mc, trial_seed) -> list(X, y,
#'                      true_support). Called once per trial.
#' @param ctrl          trex_control() object (same for all trials).
#' @param tFDR          Target FDR level.
#' @param base_seed     Base RNG seed; trial mc uses base_seed + mc - 1.
#' @param num_MC        Number of Monte Carlo trials.
#' @param num_cores     Number of parallel worker cores.
#' @param track_L_T     If TRUE, also record L and T_stop per trial.
#' @param label         Progress label printed before / after the run.
#'
#' @return Named list: mean_FDR, mean_TPR, sd_FDR, sd_TPR, and (if track_L_T) mean_L, mean_T.
.run_mc_trex <- function(make_data_fn,
                         ctrl,
                         tFDR,
                         base_seed,
                         num_MC,
                         num_cores,
                         track_L_T = FALSE,
                         label     = "") {

  one_trial <- function(mc) {
    trial_seed <- base_seed + mc - 1L
    dat <- make_data_fn(mc, trial_seed)
    sel <- TRexSelector$new(
      X       = dat$X,
      y       = dat$y,
      tFDR    = tFDR,
      seed    = -1L,  # hardware entropy per trial (indepedent random experiments)
      verbose = FALSE,
      control = ctrl
    )
    sel$select()
    res <- c(
      FDP = TRexSelectorNeo::compute_fdp(sel$selected_indices, dat$true_support),
      TPP = TRexSelectorNeo::compute_tpp(sel$selected_indices, dat$true_support)
    )
    if (track_L_T) res <- c(res, L = sel$L, T = sel$T_stop)
    res
  }

  if (nchar(label) > 0L)
    cat(sprintf("  %s \u2014 Running %d MC trials ...\n", label, num_MC))

  if (.Platform$OS.type == "unix") {
    # NOTE: mclappy on Unix forks a new R process for each trial
    # mc.preschedule = FALSE: each trial runs in its own forked child process.
    # When the child exits the OS reclaims all memory (R heap + C++ heap)
    # unconditionally — no dependency on R's GC or finalizer timing.
    #
    # This mirrors the C++ behaviour where TRexSelector is stack-allocated and
    # destroyed deterministically at end of the omp parallel for loop body.
    #
    # With mc.preschedule = TRUE (the default), long-lived workers accumulate
    # C++ stored_dummies_ across trials because R6's internal self-reference
    # cycle prevents gc() from reliably firing the EXTPTR finalizer in time.
    res_list <- parallel::mclapply(
      X               = seq_len(num_MC),
      FUN             = one_trial,
      mc.cores        = num_cores,
      mc.preschedule  = FALSE
    )
  } else {
    # NOTE: parLapply on Windows spawns a new R process for each trial
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(TRexSelectorNeo))
    parallel::clusterExport(
      cl,
      varlist = c("base_seed", "make_data_fn", "ctrl", "tFDR", "track_L_T"),
      envir   = environment()
    )
    res_list <- parallel::parLapply(cl, seq_len(num_MC), one_trial)
  }

  res_mat <- do.call(rbind, res_list)
  fdps <- res_mat[, "FDP"]
  tpps <- res_mat[, "TPP"]

  out <- list(
    mean_FDR = mean(fdps),
    mean_TPR = mean(tpps),
    sd_FDR   = sd(fdps),
    sd_TPR   = sd(tpps)
  )
  if (track_L_T) {
    out$mean_L <- mean(res_mat[, "L"])
    out$mean_T <- mean(res_mat[, "T"])
  }

  if (nchar(label) > 0L)
    cat(sprintf("  %s \u2014 done. mean_TPR=%.3f  mean_FDR=%.3f\n\n",
                label, out$mean_TPR, out$mean_FDR))
  out
}


# ==============================================================================
# .print_solver_table
# ==============================================================================

#' Print an aligned solver x SNR results table to the console.
#'
#' Columns are SNR values; rows are (solver, metric) pairs.
#'
#' @param solver_names   Character vector of solver names (row order).
#' @param snr_values     Numeric vector of SNR grid values (column order).
#' @param fdr_mat        Matrix (n_solvers x n_snr) of mean FDR values.
#' @param tpr_mat        Matrix (n_solvers x n_snr) of mean TPR values.
#' @param avg_L_mat      Matrix or NULL; if non-NULL, print avg L rows.
#' @param avg_T_mat      Matrix or NULL; if non-NULL, print avg T rows.
#' @param num_MC         Number of MC trials (printed in header).
.print_solver_table <- function(solver_names,
                                snr_values,
                                fdr_mat,
                                tpr_mat,
                                avg_L_mat = NULL,
                                avg_T_mat = NULL,
                                num_MC    = NULL) {

  # Solver column grows to fit the longest name; 15 keeps the classic layout
  # for demos whose names all fit (mirrors trex_sim_utils.hpp).
  sw <- max(15L, max(nchar(solver_names)) + 2L)
  mw <- 8L    # metric column width
  nw <- 5L    # SNR label column width
  cw <- 10L   # data column width

  n_snr <- length(snr_values)

  if (!is.null(num_MC)) {
    cat(sprintf(
      "\n%s\n=== T-Rex Results (averaged over %d Monte Carlo runs) ===\n%s\n\n",
      strrep("=", 70L), num_MC, strrep("=", 70L)
    ))
  }

  # Header
  hdr <- sprintf("%-*s%-*s%*s", sw, "Solver", mw, "Metric", nw, "SNR")
  for (snr in snr_values) {
    hdr <- paste0(hdr, sprintf("%*.1f", cw, snr))
  }

  cat(hdr, "\n")
  cat(strrep("-", sw + mw + nw + cw * n_snr), "\n")

  print_row <- function(sname, metric, vals, first_row) {

    row <- sprintf("%-*s%-*s%*s",
                   sw, if (first_row) sname else "",
                   mw, metric,
                   nw, "")

    for (v in vals) {
      row <- paste0(row, sprintf("%*.4f", cw, v))
    }
    cat(row, "\n")
  }

  for (i in seq_along(solver_names)) {
    nm <- solver_names[[i]]
    print_row(nm, "FDR", fdr_mat[i, ], TRUE)
    print_row(nm, "TPR", tpr_mat[i, ], FALSE)
    if (!is.null(avg_L_mat)) print_row(nm, "Avg L", avg_L_mat[i, ], FALSE)
    if (!is.null(avg_T_mat)) print_row(nm, "Avg T", avg_T_mat[i, ], FALSE)
    cat("\n")
  }
}


# ==============================================================================
# .save_mc_csv
# ==============================================================================

#' Write a tidy long-format CSV (solver, metric, snr, value).
#'
#' @param filepath     Full path to the output CSV file.
#' @param solver_names Character vector of solver names.
#' @param snr_values   Numeric vector of SNR grid values.
#' @param fdr_mat      Matrix (n_solvers x n_snr) of mean FDR values.
#' @param tpr_mat      Matrix (n_solvers x n_snr) of mean TPR values.
#' @param avg_L_mat    Matrix or NULL.
#' @param avg_T_mat    Matrix or NULL.
.save_mc_csv <- function(filepath,
                         solver_names,
                         snr_values,
                         fdr_mat,
                         tpr_mat,
                         avg_L_mat = NULL,
                         avg_T_mat = NULL) {

  rows <- list()
  for (i in seq_along(solver_names)) {
    nm <- solver_names[[i]]
    for (j in seq_along(snr_values)) {

      snr_str <- sprintf("%.6f", snr_values[[j]])
      rows <- c(rows,
        list(data.frame(solver = nm, metric = "FDR", snr = snr_str,
                        value = fdr_mat[i, j], stringsAsFactors = FALSE)),
        list(data.frame(solver = nm, metric = "TPR", snr = snr_str,
                        value = tpr_mat[i, j], stringsAsFactors = FALSE))
      )

      if (!is.null(avg_L_mat))
        rows <- c(rows, list(data.frame(solver = nm, metric = "AvgL",
                                        snr = snr_str, value = avg_L_mat[i, j],
                                        stringsAsFactors = FALSE)))
      if (!is.null(avg_T_mat))
        rows <- c(rows, list(data.frame(solver = nm, metric = "AvgT",
                                        snr = snr_str, value = avg_T_mat[i, j],
                                        stringsAsFactors = FALSE)))
    }
  }

  out_df <- do.call(rbind, rows)
  utils::write.csv(out_df, filepath, row.names = FALSE, quote = FALSE)
  cat(sprintf("[Info] CSV results saved to: %s\n", filepath))
}


# ==============================================================================
# .save_and_print_mc
# ==============================================================================

#' Print results table (console + .txt file) and write tidy CSV.
#'
#' @param out_dir      Output directory (must already exist).
#' @param file_stem    Output file base name (no directory, no extension).
#' @param num_MC       Number of MC trials.
#' @param solver_names Character vector of solver names.
#' @param snr_values   Numeric vector of SNR grid values.
#' @param fdr_mat      Matrix (n_solvers x n_snr) of mean FDR values.
#' @param tpr_mat      Matrix (n_solvers x n_snr) of mean TPR values.
#' @param avg_L_mat    Matrix or NULL.
#' @param avg_T_mat    Matrix or NULL.
.save_and_print_mc <- function(out_dir,
                               file_stem,
                               num_MC,
                               solver_names,
                               snr_values,
                               fdr_mat,
                               tpr_mat,
                               avg_L_mat = NULL,
                               avg_T_mat = NULL) {

  # Write to both console and text file via connection trick
  txt_path <- file.path(out_dir, paste0(file_stem, ".txt"))
  con      <- file(txt_path, open = "wt")
  sink(con, split = TRUE)
  on.exit({
    sink()
    close(con)
  }, add = TRUE)

  .print_solver_table(solver_names, snr_values, fdr_mat, tpr_mat,
                      avg_L_mat, avg_T_mat, num_MC)

  sink()
  close(con)
  on.exit(NULL)
  cat(sprintf("[Info] Results saved to:     %s\n", txt_path))

  .save_mc_csv(
    filepath     = file.path(out_dir, paste0(file_stem, ".csv")),
    solver_names = solver_names,
    snr_values   = snr_values,
    fdr_mat      = fdr_mat,
    tpr_mat      = tpr_mat,
    avg_L_mat    = avg_L_mat,
    avg_T_mat    = avg_T_mat
  )
}

# ==============================================================================
