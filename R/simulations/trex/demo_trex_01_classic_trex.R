# ==============================================================================
# demo_trex_01_classic_trex.R
# ==============================================================================
#
# Mirror of cpp/trex_selector_methods/trex/demo_trex_01_classic_trex_inmem.cpp
#
#  Part 1: Single-run basic T-Rex demo.
#            - Low- and high-dimensional settings.
#            - Fixed true support; TLARS solver.
#            - Reports phi_prime, phi_mat, fdp_hat_mat, r_mat, voting_grid.
#
#  Part 2: MC simulation with fixed true support (solver comparison).
#            - 14 solvers; SNR sweep {0.1, 0.5, 1.0, 2.0, 5.0}.
#            - Fixed support drawn once with seed 24.
#            - Reports averaged FDR / TPR per solver x SNR.
#
#  Part 3: MC simulation with variable true support (STANDARD SIMULATION).
#            - 14 solvers; SNR sweep {0.1, 0.2, ..., 2.0, 5.0} (21 values).
#            - Support and (optionally) coefficients redrawn each trial.
#            - Reports averaged FDR / TPR / Avg L / Avg T per solver x SNR.
#
# ==============================================================================

library(TRexSelector)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) "."
)
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_gauss_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))
source(file.path(this_dir_, "trex_sim_common.R"))


# ==============================================================================
# Parameters
# ==============================================================================

OUT_DIR <- file.path(this_dir_, "simulation_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NUM_MC <- 500L

# Run flags (mirrors C++ main() — set to TRUE to execute each part)
run_part_1 <- FALSE
run_part_2 <- FALSE
run_part_3 <- TRUE  # STANDARD SIMULATION (high-dim, rnd_coef = FALSE)

# Fixed true support for Parts 1 and 2 (0-based in C++, converted to 1-based)
# C++: {27, 149, 43, 128, 42, 4}  →  R (1-based): c(28, 150, 44, 129, 43, 5)
TRUE_SUPPORT_P1 <- c(28L, 150L, 44L, 129L, 43L, 5L)


# ==============================================================================
# Shared helpers
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


# ==============================================================================
# Part 1: Single-run basic demo
# ==============================================================================
# Mirrors demo_TRexSelector() in demo_trex_01_classic_trex_inmem.cpp.
# C++ main(): if (false) demo_TRexSelector(high_dim=true, rnd_coef=false)

trx_01_part_1 <- function() {

  # Parameters matching C++ demo_TRexSelector()
  tFDR <- 0.1
  snr  <- 1.0
  # C++ true_coefs = all 1  →  amplitude = 1 (default)

  ctrl <- trex_control(
    method                   = "TLARS",
    K                        = 20L,
    max_dummy_multiplier     = 10L,
    use_max_T_stop           = TRUE,
    lloop_strategy           = "HCONCAT",
    tloop_stagnation_stop    = FALSE,
    tloop_max_stagnant_steps = 5L,
    use_openmp               = FALSE
  )

  run_single <- function(high_dim) {
    n <- if (high_dim) 150L else 5000L
    p <- if (high_dim) 300L else 1000L
    dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"

    cat("\n", strrep("=", 70L), "\n", sep = "")
    cat("Part 1: Single-run basic T-Rex demo  |  ", dim_label, "\n", sep = "")
    cat(sprintf("n = %d,  p = %d,  tFDR = %.2f,  SNR = %.1f\n\n", n, p, tFDR, snr))

    cat("Generating synthetic data ...\n")
    dat <- dgp_gauss_snr(n, p, TRUE_SUPPORT_P1, amplitude = 1, snr = snr, seed = 58L)

    cat("Running T-Rex Selector (TLARS) ...\n\n")
    selector <- TRexSelector$new(
      X       = dat$X,
      y       = dat$y,
      tFDR    = tFDR,
      seed    = 58L,
      verbose = TRUE,
      control = ctrl
    )
    selector$select()

    stem <- sprintf("d01_trex_basic_n%d_p%d", n, p)
    .print_single_run(
      sprintf("Part 1 \u2014 i.i.d. Gaussian  [%s, TLARS]", dim_label),
      dat, selector, tFDR
    )
    invisible(selector)
  }

  run_single(high_dim = FALSE)
  run_single(high_dim = TRUE)
}


# ==============================================================================
# Part 2: MC simulation — fixed true support, solver comparison
# ==============================================================================
# Mirrors demo_TRexSelector_MonteCarlo() in demo_trex_01_classic_trex_inmem.cpp.

trx_01_part_2 <- function(num_MC = NUM_MC, high_dim = TRUE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR <- 0.1
  snr_values   <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  stagnant_sw  <- 7L

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 2: MC Simulation — Fixed True Support  |  ", dim_label, "\n", sep = "")
  cat(sprintf("n = %d,  p = %d,  s = %d,  num_MC = %d\n", n, p, s, num_MC))

  # Base control (per-solver params added in .make_trex_ctrl)
  base_ctrl <- list(
    K                        = 20L,
    max_dummy_multiplier     = 10L,
    use_max_T_stop           = TRUE,
    lloop_strategy           = "HCONCAT",
    tloop_stagnation_stop    = TRUE,
    tloop_max_stagnant_steps = stagnant_sw,
    use_openmp               = FALSE,
    opt_threshold            = 0.75,
    use_memory_mapping       = FALSE
  )

  # Generate the fixed true support once (shared across all MC runs)
  # Mirrors: std::mt19937 rng(24); udist(0, p-1); cardinality=10
  true_support <- make_support_random(s, p, seed = 24L)
  cat(sprintf("True support (cardinality %d): {%s}\n\n",
              s, paste(true_support, collapse = ", ")))

  # True coefficients: all 1 (rnd_coef = FALSE, matching C++ default)
  true_coefs <- rep(1.0, s)

  n_sol <- length(SOLVERS_DEFAULT)
  n_snr <- length(snr_values)
  fdr_mat <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  tpr_mat <- matrix(0.0, nrow = n_sol, ncol = n_snr)
  rownames(fdr_mat) <- rownames(tpr_mat) <- vapply(SOLVERS_DEFAULT, `[[`, "", "name")
  colnames(fdr_mat) <- colnames(tpr_mat) <- as.character(snr_values)

  for (si in seq_along(SOLVERS_DEFAULT)) {
    solver <- SOLVERS_DEFAULT[[si]]
    ctrl   <- .make_trex_ctrl(solver, base_ctrl)

    cat(sprintf("Solver: %s\n", solver$name))

    for (snr_idx in seq_along(snr_values)) {
      snr        <- snr_values[[snr_idx]]
      base_seed  <- 24L + (snr_idx - 1L) * 1000L

      make_data <- function(mc, trial_seed) {
        dgp_gauss_snr(n, p, true_support, coefs = true_coefs, snr = snr,
                      seed = trial_seed)
      }

      lbl    <- sprintf("SNR=%.2f [%s]", snr, solver$name)
      result <- .run_mc_trex(make_data, ctrl, tFDR, base_seed, num_MC,
                             num_cores, track_L_T = FALSE, label = lbl)

      fdr_mat[si, snr_idx] <- result$mean_FDR
      tpr_mat[si, snr_idx] <- result$mean_TPR
    }
    cat("\n")
  }

  stem <- sprintf("demo_trex_01_p02_trex_results_n%d_p%d_stagnation_window_%d",
                  n, p, stagnant_sw)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     rownames(fdr_mat), snr_values, fdr_mat, tpr_mat)
  cat("\n")
}


# ==============================================================================
# Part 3: MC simulation — variable true support (STANDARD SIMULATION)
# ==============================================================================
# Mirrors demo_TRexSelector_varMonteCarlo() in demo_trex_01_classic_trex_inmem.cpp.
# C++ main(): if (true) demo_TRexSelector_varMonteCarlo(500, high_dim=true, rnd_coef=false)

trx_01_part_3 <- function(num_MC = NUM_MC, high_dim = TRUE, rnd_coef = FALSE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR <- 0.1
  # SNR = {0.1, 0.2, ..., 2.0, 5.0}  (21 values — mirrors C++ std::iota trick)
  snr_values   <- c(seq(0.1, 2.0, by = 0.1), 5.0)
  stagnant_sw  <- 3L

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 3: MC Simulation — Variable Support  |  ", dim_label, "\n", sep = "")
  cat(sprintf("n = %d,  p = %d,  s = %d,  num_MC = %d\n\n", n, p, s, num_MC))

  # Base control
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

  for (si in seq_along(SOLVERS_DEFAULT)) {
    solver <- SOLVERS_DEFAULT[[si]]
    ctrl   <- .make_trex_ctrl(solver, base_ctrl)

    cat(sprintf("%s\nSolver: %s\n%s\n\n",
                strrep("=", 50L), solver$name, strrep("=", 50L)))

    for (snr_idx in seq_along(snr_values)) {
      snr       <- snr_values[[snr_idx]]
      base_seed <- 24L + (snr_idx - 1L) * 1000L

      make_data <- function(mc, trial_seed) {
        # Per-trial support: mirrors C++ rng_sup(seed + 500000u)
        trial_support <- make_support_random(s, p, seed = trial_seed)
        # Per-trial coefficients: all 1 (rnd_coef = FALSE)
        trial_coefs <- if (rnd_coef) {
          set.seed(trial_seed + 600000L)
          stats::rnorm(s)
        } else {
          rep(1.0, s)
        }
        dgp_gauss_snr(n, p, trial_support, coefs = trial_coefs, snr = snr,
                      seed = trial_seed)
      }

      lbl    <- sprintf("SNR=%.2f [%s]", snr, solver$name)
      result <- .run_mc_trex(make_data, ctrl, tFDR, base_seed, num_MC,
                             num_cores, track_L_T = TRUE, label = lbl)

      fdr_mat[si, snr_idx]   <- result$mean_FDR
      tpr_mat[si, snr_idx]   <- result$mean_TPR
      avg_L_mat[si, snr_idx] <- result$mean_L
      avg_T_mat[si, snr_idx] <- result$mean_T
    }
    cat("\n")
  }

  stem <- sprintf("demo_trex_01_p03_trex_results_n%d_p%d_stagnation_window_%d",
                  n, p, stagnant_sw)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     sol_names, snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_1) trx_01_part_1()
if (run_part_2) trx_01_part_2(NUM_MC, high_dim = TRUE)
if (run_part_3) trx_01_part_3(NUM_MC, high_dim = TRUE, rnd_coef = FALSE)

cat("\ndemo_trex_01_classic_trex.R complete.\n")
# ==============================================================================
