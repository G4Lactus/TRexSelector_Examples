# ==============================================================================
# trex_helpers.R
# ==============================================================================
# Shared helper functions for TRex demos and simulations.
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
  tpp          <- TRexSelector::compute_tpp(selected, true_support)
  fdp          <- TRexSelector::compute_fdp(selected, true_support)

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
