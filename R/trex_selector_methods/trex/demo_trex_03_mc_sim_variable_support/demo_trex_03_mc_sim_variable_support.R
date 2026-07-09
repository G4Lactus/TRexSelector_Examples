# ==============================================================================
# demo_trex_03_mc_sim_variable_support.R
# ==============================================================================
#
# Classic T-Rex selector demo, where in part 3 a Monte Carlo (MC) simulation
# with variable true support and sweep over all 14 tsolvers is performed.
#
# The file mirrors:
# cpp/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support.cpp
#
# Demo content:
# ----------------------------------------------------------------
#
#  - 14 solvers; SNR sweep {0.1, 0.2, ..., 2.0, 5.0} (21 values).
#  - Support and (optionally) coefficients redrawn each trial.
#  - Reports averaged FDR / TPR / Avg L / Avg T per solver x SNR.
#
# ==============================================================================

library(TRexSelectorNeo)
library(parallel)

# ==============================================================================
# Setup
# ==============================================================================
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

source(file.path(this_dir_, "..", "..", "support_generators.R"))
source(file.path(this_dir_, "..", "..", "simulation_utils.R"))
source(file.path(this_dir_, "..", "trex_sim_utils.R"))


# ==============================================================================
# Global Simulation Parameters
# ==============================================================================

OUT_DIR <- file.path(this_dir_, "simulation_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

num_cores <- 6L
NUM_MC <- 200L
run_trex_03 <- TRUE


# ==============================================================================
# MC simulation — variable true support (STANDARD SIMULATION)
# ==============================================================================

trx_03_mc_sim_variable_support <- function(num_MC = NUM_MC, high_dim = TRUE, rnd_coef = FALSE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L

  s <- 10L
  tFDR <- 0.1
  snr_values   <- c(seq(0.1, 2.0, by = 0.1), 5.0)
  stagnant_sw  <- 5L

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("MC Simulation \u2014 Variable Support  |  ", dim_label, "\n", sep = "")
  cat(sprintf("n = %d,  p = %d,  s = %d,  num_MC = %d\n\n", n, p, s, num_MC))

  base_ctrl <- list(
    K                        = 20L,
    max_dummy_multiplier     = 10L,
    use_max_T_stop           = TRUE,
    lloop_strategy           = "HCONCAT",
    tloop_stagnation_stop    = TRUE,
    tloop_max_stagnant_steps = stagnant_sw,
    use_openmp               = FALSE
  )

  n_sol <- length(SOLVERS_DEFAULT)
  n_snr <- length(snr_values)
  fdr_mat   <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  tpr_mat   <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  avg_L_mat <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  avg_T_mat <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  sol_names <- vapply(SOLVERS_DEFAULT, `[[`, "", "name")
  rownames(fdr_mat) <- rownames(tpr_mat) <-
    rownames(avg_L_mat) <- rownames(avg_T_mat) <- sol_names
  colnames(fdr_mat) <- colnames(tpr_mat) <-
    colnames(avg_L_mat) <- colnames(avg_T_mat) <- as.character(snr_values)

  # Simulation loop: iterate over solvers and SNR values
  for (si in seq_along(SOLVERS_DEFAULT)) {
    solver <- SOLVERS_DEFAULT[[si]]
    ctrl   <- .make_trex_ctrl(solver, base_ctrl)
    cat(sprintf("%s\nSolver: %s\n%s\n\n", strrep("=", 50L), solver$name, strrep("=", 50L)))
    for (snr_idx in seq_along(snr_values)) {
      snr       <- snr_values[[snr_idx]]
      base_seed <- 24L + (snr_idx - 1L) * 1000L

      make_data <- local({
        current_snr <- snr
        function(mc, trial_seed) {
          trial_support <- make_support_random(s, p, seed = trial_seed)
          trial_coefs <- if (rnd_coef) {
            # Use an isolated seed stream
            set.seed(trial_seed + 600000L)
            stats::rnorm(s)
          } else {
            rep(1.0, s)
          }
          dgp_gauss_snr(n, p, trial_support, coefs = trial_coefs,
                        snr = current_snr, seed = trial_seed)
        }
      })

      lbl    <- sprintf("SNR=%.2f [%s]", snr, solver$name)
      result <- .run_mc_trex(make_data, ctrl, tFDR, base_seed, num_MC, num_cores,
                             track_L_T = TRUE, label = lbl)

      fdr_mat[si, snr_idx]   <- result$mean_FDR
      tpr_mat[si, snr_idx]   <- result$mean_TPR
      avg_L_mat[si, snr_idx] <- result$mean_L
      avg_T_mat[si, snr_idx] <- result$mean_T
    }
    cat("\n")
  }

  # Save results
  stem <- sprintf("demo_trex_03_mc_sim_variable_support_results_n%d_p%d_stagnation_window_%d",
                  n, p, stagnant_sw)
  .save_and_print_mc(OUT_DIR, stem, num_MC, sol_names, snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Main
# -------------
# Run Part 3
# ==============================================================================

if (run_trex_03) trx_03_mc_sim_variable_support(NUM_MC, high_dim = TRUE, rnd_coef = FALSE)
cat("\ndemo_trex_03_mc_sim_variable_support.R complete.\n")
