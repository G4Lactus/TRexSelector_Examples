# ==============================================================================
# trex_sim_common.R
# ==============================================================================
#
# Shared infrastructure for classic T-Rex MC simulations.
# Sourced by demo_trex_0{1,2,3}_*.R — do not run standalone.
#
# Mirrors C++ trex_sim_common.hpp.
#
# Contents:
#   SOLVERS_DEFAULT        — list of 14 solver descriptors (mirrors
#                            make_default_solvers_to_test() in C++).
#   .make_trex_ctrl()      — build trex_control() from a solver descriptor and
#                            a named base-params list.
#   .run_mc_trex()         — parallel MC runner using the TRexSelector OOP API.
#                            Tracks FDP, TPP and optionally L, T_stop.
#   .print_solver_table()  — print aligned solver x SNR table to console.
#   .save_mc_csv()         — write tidy long-format CSV to disk.
#   .save_and_print_mc()   — run .print_solver_table() + .save_mc_csv() + .txt.
#
# ==============================================================================


# ==============================================================================
# Default solver list
# (mirrors make_default_solvers_to_test() in demo_trex_01_classic_trex_inmem.cpp)
# ==============================================================================

SOLVERS_DEFAULT <- list(
  list(method = "TLARS",      name = "TLARS",        lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TLASSO",     name = "TLASSO",       lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TENET",      name = "TENET",        lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TSTAGEWISE", name = "TSTAGEWISE",   lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TSTEPWISE",  name = "TSTEPWISE",    lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TOMP",       name = "TOMP",         lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TGP",        name = "TGP",          lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TACGP",      name = "TACGP",        lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TMP",        name = "TMP",          lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TAFS",       name = "TAFS_rho_0.3", lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TAFS",       name = "TAFS_rho_1.0", lambda2 = 0.1, rho_afs = 1.0,  ncgmp_variant = 0L),
  list(method = "TNCGMP",     name = "TNCGMP_v1",    lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 1L),
  list(method = "TNCGMP",     name = "TNCGMP_v0",    lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L),
  list(method = "TOOLS",      name = "TOOLS",        lambda2 = 0.1, rho_afs = 0.3,  ncgmp_variant = 0L)
)


# ==============================================================================
# .make_trex_ctrl
# ==============================================================================

#' Build a trex_control() object from a solver descriptor and a base-params list.
#'
#' @param solver    One element of SOLVERS_DEFAULT (a named list with method,
#'                  lambda2, rho_afs, ncgmp_variant).
#' @param base      Named list of shared control parameters forwarded verbatim
#'                  to trex_control().  Must NOT include method, lambda2,
#'                  rho_afs or ncgmp_variant (those come from `solver`).
#'
#' @return A trex_control() object.
.make_trex_ctrl <- function(solver, base) {
  do.call(
    trex_control,
    c(list(method        = solver$method,
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
#'                        true_support).  Called once per trial.
#' @param ctrl          trex_control() object (same for all trials).
#' @param tFDR          Target FDR level.
#' @param base_seed     Base RNG seed; trial mc uses base_seed + mc - 1.
#' @param num_MC        Number of Monte Carlo trials.
#' @param num_cores     Number of parallel worker cores.
#' @param track_L_T     If TRUE, also record L and T_stop per trial.
#' @param label         Progress label printed before / after the run.
#'
#' @return Named list: mean_FDR, mean_TPR, sd_FDR, sd_TPR, and (if
#'   track_L_T) mean_L, mean_T.
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
      seed    = trial_seed,
      verbose = FALSE,
      control = ctrl
    )
    sel$select()
    res <- c(
      FDP = TRexSelector::compute_fdp(sel$selected_indices, dat$true_support),
      TPP = TRexSelector::compute_tpp(sel$selected_indices, dat$true_support)
    )
    if (track_L_T) res <- c(res, L = sel$L, T = sel$T_stop)
    res
  }

  if (nchar(label) > 0L)
    cat(sprintf("  %s \u2014 Running %d MC trials ...\n", label, num_MC))

  if (.Platform$OS.type == "unix") {
    res_list <- parallel::mclapply(
      X        = seq_len(num_MC),
      FUN      = one_trial,
      mc.cores = num_cores
    )
  } else {
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(TRexSelector))
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

  sw <- 15L   # solver column width
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
  for (snr in snr_values)
    hdr <- paste0(hdr, sprintf("%*.1f", cw, snr))
  cat(hdr, "\n")
  cat(strrep("-", sw + mw + nw + cw * n_snr), "\n")

  print_row <- function(sname, metric, vals, first_row) {
    row <- sprintf("%-*s%-*s%*s",
                   sw, if (first_row) sname else "",
                   mw, metric,
                   nw, "")
    for (v in vals)
      row <- paste0(row, sprintf("%*.4f", cw, v))
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
  on.exit({ sink(); close(con) }, add = TRUE)

  .print_solver_table(solver_names, snr_values, fdr_mat, tpr_mat,
                      avg_L_mat, avg_T_mat, num_MC)

  sink(); close(con)
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
