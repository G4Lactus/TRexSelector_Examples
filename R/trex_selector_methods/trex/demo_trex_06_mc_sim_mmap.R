# ==============================================================================
# demo_trex_06_mmap_mc_sim.R
# ==============================================================================
#
# Classical T-Rex selector demo.
# Part 6 is a Monte Carlo (MC) simulation with memory-mapped X and
# use_memory_mapping = TRUE to verify the D-mmap + solver-serialization.
#
# The file mirrors:
# cpp/trex_selector_methods/trex/demo_trex_06_mmap_mc_sim.cpp
#
# Demo content:
# ----------------------------------------------------------------
#
#  Part A: MC — In-memory X + use_memory_mapping = TRUE.
#            Mirrors C++ Demo A.
#            Activates D-mmap + solver serialization inside TRexSelector.
#            Purpose: verify the D-mmap + solver-serialization pipeline
#            produces stable, reproducible FDR/TPR statistics over many runs.
#
#  Part B: MC — Memory-mapped X + use_memory_mapping = TRUE.
#            Mirrors C++ Demo B.
#            Each MC trial creates a temporary mmap-backed X/y file and
#            removes it via on.exit() (exception-safe RAII equivalent).
#
# Both parts use:
#   - only the TLARS solver is demonstrated here
#   - Fixed true support (drawn once with seed 24)
#   - SNR sweep: {0.1, 0.5, 1.0, 2.0, 5.0}
#   - tloop_max_stagnant_steps = 7,  opt_threshold = 0.75
#
# ==============================================================================

library(TRexSelector)
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
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_gauss_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))
source(file.path(this_dir_, "trex_sim_common.R"))


# ==============================================================================
# Global Simulation Parameters
# ==============================================================================
OUT_DIR <- file.path(this_dir_, "simulation_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

num_cores <- 6L
NUM_MC <- 500L

# Run flags
run_part_a <- TRUE
run_part_b <- TRUE


# ==============================================================================
# Part A: MC — In-memory X + use_memory_mapping = TRUE
# ==============================================================================
# Mirrors demo_TRexSelector_d_mmap_MonteCarlo() in C++.

trex_06_part_a <- function(num_MC = NUM_MC, high_dim = TRUE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR       <- 0.1
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)

  .print_demo_header(
    "Part C: MC \u2014 In-Memory X + use_memory_mapping = TRUE",
    "Verifies D-mmap + solver-serialization pipeline over many runs.",
    high_dim, n, p, num_MC
  )

  ctrl <- make_mmap_trex_ctrl(
    solver                   = "TLARS",
    tloop_max_stagnant_steps = 7L,
    use_openmp               = FALSE,
    opt_threshold            = 0.75
  )

  # Fixed true support shared across all MC runs (mirrors C++ rng(24))
  true_support <- make_support_random(s, p, seed = 24L)
  true_coefs   <- rep(1.0, s)   # rnd_coef = FALSE
  cat(sprintf("True support (cardinality %d): {%s}\n\n",
              s, paste(true_support, collapse = ", ")))

  n_snr     <- length(snr_values)
  fdr_vec   <- numeric(n_snr)
  tpr_vec   <- numeric(n_snr)
  avg_L_vec <- numeric(n_snr)
  avg_T_vec <- numeric(n_snr)

  for (snr_idx in seq_along(snr_values)) {
    snr       <- snr_values[[snr_idx]]
    base_seed <- 24L + (snr_idx - 1L) * 1000L

    make_data <- function(mc, trial_seed) {
      dgp_gauss_snr(n, p, true_support, coefs = true_coefs, snr = snr,
                    seed = trial_seed)
    }

    lbl    <- sprintf("SNR=%.2f [TLARS/mmap-D]", snr)
    result <- .run_mc_trex(make_data, ctrl, tFDR, base_seed, num_MC,
                           num_cores, track_L_T = TRUE, label = lbl)

    fdr_vec[[snr_idx]]   <- result$mean_FDR
    tpr_vec[[snr_idx]]   <- result$mean_TPR
    avg_L_vec[[snr_idx]] <- result$mean_L
    avg_T_vec[[snr_idx]] <- result$mean_T
  }

  # Build 1-row result matrices for reuse of .save_and_print_mc()
  fdr_mat   <- matrix(fdr_vec,   nrow = 1L, dimnames = list("TLARS", NULL))
  tpr_mat   <- matrix(tpr_vec,   nrow = 1L, dimnames = list("TLARS", NULL))
  avg_L_mat <- matrix(avg_L_vec, nrow = 1L, dimnames = list("TLARS", NULL))
  avg_T_mat <- matrix(avg_T_vec, nrow = 1L, dimnames = list("TLARS", NULL))

  stem <- sprintf("d03_trex_mmap_demo_c_n%d_p%d_sw%d", n, p,
                  ctrl$tloop_max_stagnant_steps)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     "TLARS", snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Part B: MC — Memory-mapped X + use_memory_mapping = TRUE
# ==============================================================================
# Mirrors demo_TRexSelector_full_mmap_MonteCarlo() in C++.
#
# Each MC trial writes X to a per-trial temporary mmap file and removes it
# via on.exit() inside the factory lambda — mirroring the C++ RAII MmapFileGuard
# declared before SyntheticDataMapped (LIFO cleanup, exception-safe).
#
# NOTE: mclapply forks child processes, so each trial runs in its own address
# space; no cross-trial file conflicts are possible when using tempfile().

trex_06_part_b <- function(num_MC = NUM_MC, high_dim = TRUE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR       <- 0.1
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)

  .print_demo_header(
    "Part D: MC \u2014 Memory-Mapped X + use_memory_mapping = TRUE",
    "Full mmap pipeline: X mmap + D mmap + solver serialization.",
    high_dim, n, p, num_MC
  )

  ctrl <- make_mmap_trex_ctrl(
    solver                   = "TLARS",
    tloop_max_stagnant_steps = 7L,
    use_openmp               = FALSE,
    opt_threshold            = 0.75
  )

  # Fixed true support shared across all MC runs
  true_support <- make_support_random(s, p, seed = 24L)
  true_coefs   <- rep(1.0, s)
  cat(sprintf("True support (cardinality %d): {%s}\n\n",
              s, paste(true_support, collapse = ", ")))

  n_snr     <- length(snr_values)
  fdr_vec   <- numeric(n_snr)
  tpr_vec   <- numeric(n_snr)
  avg_L_vec <- numeric(n_snr)
  avg_T_vec <- numeric(n_snr)

  for (snr_idx in seq_along(snr_values)) {
    snr       <- snr_values[[snr_idx]]
    base_seed <- 24L + (snr_idx - 1L) * 1000L

    # Factory: generates X in memory, writes to mmap file, creates
    # MemoryMappedMatrix, runs selector, removes file on exit.
    # Mirrors C++ make_data_mmap lambda with MmapFileGuard + SyntheticDataMapped.
    make_data_mmap <- function(mc, trial_seed) {
      dat    <- dgp_gauss_snr(n, p, true_support, coefs = true_coefs, snr = snr,
                              seed = trial_seed)

      # execute_with_temp_mmap safely handles the temporary file lifecycle
      X_mmap <- execute_with_temp_mmap(dat$X, function(x) x)

      # Return a list where X is the mmap object; y and true_support are copies
      list(X = X_mmap, y = dat$y, true_support = dat$true_support)
    }

    lbl    <- sprintf("SNR=%.2f [TLARS/mmap-XD]", snr)
    result <- .run_mc_trex(make_data_mmap, ctrl, tFDR, base_seed, num_MC,
                           num_cores, track_L_T = TRUE, label = lbl)

    fdr_vec[[snr_idx]]   <- result$mean_FDR
    tpr_vec[[snr_idx]]   <- result$mean_TPR
    avg_L_vec[[snr_idx]] <- result$mean_L
    avg_T_vec[[snr_idx]] <- result$mean_T
  }

  fdr_mat   <- matrix(fdr_vec,   nrow = 1L, dimnames = list("TLARS", NULL))
  tpr_mat   <- matrix(tpr_vec,   nrow = 1L, dimnames = list("TLARS", NULL))
  avg_L_mat <- matrix(avg_L_vec, nrow = 1L, dimnames = list("TLARS", NULL))
  avg_T_mat <- matrix(avg_T_vec, nrow = 1L, dimnames = list("TLARS", NULL))

  stem <- sprintf("d03_trex_mmap_demo_d_n%d_p%d_sw%d", n, p,
                  ctrl$tloop_max_stagnant_steps)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     "TLARS", snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Main
# -------------
# Run parts C and D
# ==============================================================================
if (run_part_c) trex_06_part_a(NUM_MC, high_dim = TRUE)
if (run_part_d) trex_06_part_b(NUM_MC, high_dim = TRUE)
cat("\ndemo_trex_06_mmap_mc_sim.R complete.\n")
