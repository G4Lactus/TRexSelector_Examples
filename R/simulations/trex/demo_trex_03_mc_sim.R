# ==============================================================================
# demo_trex_03_mc_sim.R
# ==============================================================================
#
# Mirror of cpp/trex_selector_methods/trex/demo_trex_03_classic_trex_mmap_mc_sim.cpp
#
# Contains two MC simulation scenarios that mirror demo_trex_02_classic_trex.R:
#
#  Part C: MC — In-memory X + use_memory_mapping = TRUE.
#            Mirrors C++ Demo C.
#            Activates D-mmap + solver serialization inside TRexSelector.
#            Purpose: verify the D-mmap + solver-serialization pipeline
#            produces stable, reproducible FDR/TPR statistics over many runs.
#
#  Part D: MC — Memory-mapped X + use_memory_mapping = TRUE.
#            Mirrors C++ Demo D.
#            Each MC trial creates a temporary mmap-backed X/y file and
#            removes it via on.exit() (exception-safe RAII equivalent).
#
# Both parts use:
#   - TLARS solver only
#   - Fixed true support (drawn once with seed 24)
#   - SNR sweep: {0.1, 0.5, 1.0, 2.0, 5.0}
#   - tloop_max_stagnant_steps = 7,  opt_threshold = 0.75
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

# Run flags
run_part_c <- TRUE
run_part_d <- TRUE

# Shared TRex control (mirrors make_mmap_trex_ctrl() in C++, same as demo_02)
# use_memory_mapping = TRUE enables D-mmap + solver serialization.
MMAP_CTRL_BASE <- list(
  method                   = "TLARS",
  K                        = 20L,
  max_dummy_multiplier     = 10L,
  use_max_T_stop           = TRUE,
  lloop_strategy           = "HCONCAT",
  tloop_stagnation_stop    = TRUE,
  tloop_max_stagnant_steps = 7L,
  use_openmp               = FALSE,
  opt_threshold            = 0.75,
  use_memory_mapping       = TRUE
)


# ==============================================================================
# Part C: MC — In-memory X + use_memory_mapping = TRUE
# ==============================================================================
# Mirrors demo_TRexSelector_d_mmap_MonteCarlo() in C++.

trx_03_part_c <- function(num_MC = NUM_MC, high_dim = TRUE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR       <- 0.1
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part C: MC \u2014 In-Memory X + use_memory_mapping = TRUE\n")
  cat("  Verifies D-mmap + solver-serialization pipeline over many runs.\n")
  cat(sprintf("  %s  (n = %d, p = %d,  num_MC = %d)\n\n",
              dim_label, n, p, num_MC))

  ctrl <- do.call(trex_control, MMAP_CTRL_BASE)

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
                  MMAP_CTRL_BASE$tloop_max_stagnant_steps)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     "TLARS", snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Part D: MC — Memory-mapped X + use_memory_mapping = TRUE
# ==============================================================================
# Mirrors demo_TRexSelector_full_mmap_MonteCarlo() in C++.
#
# Each MC trial writes X to a per-trial temporary mmap file and removes it
# via on.exit() inside the factory lambda — mirroring the C++ RAII MmapFileGuard
# declared before SyntheticDataMapped (LIFO cleanup, exception-safe).
#
# NOTE: mclapply forks child processes, so each trial runs in its own address
# space; no cross-trial file conflicts are possible when using tempfile().

trx_03_part_d <- function(num_MC = NUM_MC, high_dim = TRUE) {

  n <- if (high_dim) 300L else 1000L
  p <- if (high_dim) 1000L else 300L
  s <- 10L
  tFDR       <- 0.1
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)

  dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part D: MC \u2014 Memory-Mapped X + use_memory_mapping = TRUE\n")
  cat("  Full mmap pipeline: X mmap + D mmap + solver serialization.\n")
  cat(sprintf("  %s  (n = %d, p = %d,  num_MC = %d)\n\n",
              dim_label, n, p, num_MC))

  ctrl <- do.call(trex_control, MMAP_CTRL_BASE)

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
      x_path <- tempfile(fileext = ".dat")
      on.exit(unlink(x_path), add = TRUE)   # RAII: remove file when scope exits

      dat    <- dgp_gauss_snr(n, p, true_support, coefs = true_coefs, snr = snr,
                               seed = trial_seed)
      X_mmap <- TRexSelector::convert_to_memory_mapped(dat$X, x_path)

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
                  MMAP_CTRL_BASE$tloop_max_stagnant_steps)
  .save_and_print_mc(OUT_DIR, stem, num_MC,
                     "TLARS", snr_values, fdr_mat, tpr_mat,
                     avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_c) trx_03_part_c(NUM_MC, high_dim = TRUE)
if (run_part_d) trx_03_part_d(NUM_MC, high_dim = TRUE)

cat("\ndemo_trex_03_mc_sim.R complete.\n")
# ==============================================================================
