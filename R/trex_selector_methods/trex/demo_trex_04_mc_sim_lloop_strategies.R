# ==============================================================================
# demo_trex_04_lloop_strategies.R
# ==============================================================================
#
# Classical T-Rex selector demo, where in part 4 a Monte Carlo simulation
# compares the six L-loop strategies implemented in the TRexSelector currently.
# It is demonstrated here with the TLARS solver as TRex's base selector.
#
# Mirrors:
# cpp/trex_selector_methods/trex/demo_trex_04_lloop_strategies.cpp
#
# Demo content:
# ----------------------------------------------------------------
#
#    STANDARD          — Fresh i.i.d. dummy matrix at each L-loop iteration
#                        (conservative gold standard; highest memory cost).
#    HCONCAT           — Horizontally expand dummy columns from the same base
#                        draw.
#    PERMUTATION       — Re-use the base dummy matrix via permutations.
#    PERMUTATION_DIRECT— Seed-based permutations; no base matrix in memory.
#    DIRECT            — Seed-based i.i.d. draws; no base matrix in memory.
#    SKIPL             — Skip the L-loop entirely; always uses
#                        L = max_dummy_multiplier (no adaptive L calibration).
#
#  Scenario:
#  ----------
#  High-dimensional  (n = 300, p = 1000, s = 10).
#  SNR sweep: {0.1, 0.2, ..., 2.0, 5.0}  (21 values).
#
#  Two support scenarios:
#    random support  — redrawn per trial;
#    block support   — contiguous block {1, 2, ..., s}.
#  Reports averaged FDR / TPR / Avg L / Avg T per L-strategy x SNR.
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
NUM_MC <- 25L
run_trex_04 <- TRUE

# ==============================================================================
# L-loop strategies to compare in Part 4
# ==============================================================================

# All L-loop strategy variants to compare.
# name           : row label used in output tables and CSV.
# strategy       : string passed to trex_control(lloop_strategy = ...).
# max_dummy_mult : only specified for SKIPL variants (fixed L = max_dummy_mult * p);
#                  adaptive strategies default to 10 via .make_strat_entry() below.

.make_strat_entry <- function(name, strategy, max_dummy_mult = 10L) {
  list(name = name, strategy = strategy, max_dummy_mult = max_dummy_mult)
}

L_STRATEGIES <- list(
  .make_strat_entry("STANDARD",           "STANDARD"),
  .make_strat_entry("HCONCAT",            "HCONCAT"),
  .make_strat_entry("PERMUTATION",        "PERMUTATION"),
  .make_strat_entry("PERMUTATION_DIRECT", "PERMUTATION_DIRECT"),
  .make_strat_entry("DIRECT",             "DIRECT"),
  .make_strat_entry("SKIPL_5p",  "SKIPL",  5L),
  .make_strat_entry("SKIPL_10p", "SKIPL", 10L),
  .make_strat_entry("SKIPL_20p", "SKIPL", 20L),
  .make_strat_entry("SKIPL_50p", "SKIPL", 50L)
)


# ==============================================================================
# MC simulation — variable true support, L-strategy comparison
# ==============================================================================

trx_04_lloop_strategies <- function(num_MC = NUM_MC, rnd_coef = FALSE,
                                    block_support = FALSE) {

  n    <- 300L
  p    <- 1000L
  s    <- 10L
  tFDR <- 0.10
  snr_values <- c(seq(0.1, 2.0, by = 0.1), 5.0)

  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("MC Simulation \xe2\x80\x94 L-Strategy Comparison  |  TLARS\n")
  cat("  High-dimensional (p > n)\n")
  cat(sprintf("  n = %d,  p = %d,  s = %d,  num_MC = %d  [%s support]\n\n",
              n, p, s, num_MC, if (block_support) "block" else "random"))

  n_strat <- length(L_STRATEGIES)
  n_snr   <- length(snr_values)
  fdr_mat   <- matrix(0.0, nrow = n_strat, ncol = n_snr)
  tpr_mat   <- matrix(0.0, nrow = n_strat, ncol = n_snr)
  avg_L_mat <- matrix(0.0, nrow = n_strat, ncol = n_snr)
  avg_T_mat <- matrix(0.0, nrow = n_strat, ncol = n_snr)
  strat_names <- vapply(L_STRATEGIES, `[[`, "", "name")
  rownames(fdr_mat) <- rownames(tpr_mat) <-
    rownames(avg_L_mat) <- rownames(avg_T_mat) <- strat_names
  colnames(fdr_mat) <- colnames(tpr_mat) <-
    colnames(avg_L_mat) <- colnames(avg_T_mat) <- as.character(snr_values)


  # Simulation sweep: iterate over L-loop strategies and SNR values
  for (si in seq_along(L_STRATEGIES)) {
    strat <- L_STRATEGIES[[si]]

    ctrl <- TRexSelector::trex_control(
      solver               = "TLARS",
      K                    = 20L,
      max_dummy_multiplier = strat$max_dummy_mult,
      use_max_T_stop       = TRUE,
      lloop_strategy       = strat$strategy
    )

    cat(sprintf("%s\nL-strategy: %s\n%s\n\n",
                strrep("=", 50L), strat$name, strrep("=", 50L)))

    for (snr_idx in seq_along(snr_values)) {
      snr       <- snr_values[[snr_idx]]
      base_seed <- 24L + (snr_idx - 1L) * 1000L

      make_data <- local({
        current_snr <- snr
        use_block   <- block_support
        function(mc, trial_seed) {
          trial_support <- if (use_block) seq_len(s) else
            make_support_random(s, p, seed = trial_seed)
          trial_coefs <- rep(1.0, s)
          dgp_gauss_snr(n, p, trial_support, coefs = trial_coefs,
                        snr = current_snr, seed = trial_seed)
        }
      })

      lbl    <- sprintf("SNR=%.2f [%s]", snr, strat$name)
      result <- .run_mc_trex(make_data, ctrl, tFDR, base_seed, num_MC, num_cores,
                             track_L_T = TRUE, label = lbl)

      fdr_mat[si, snr_idx]   <- result$mean_FDR
      tpr_mat[si, snr_idx]   <- result$mean_TPR
      avg_L_mat[si, snr_idx] <- result$mean_L
      avg_T_mat[si, snr_idx] <- result$mean_T
    }
    cat("\n")
  }

  stem <- sprintf("demo_trex_04_lloop_strategies_results_n%d_p%d_%s",
                  n, p, if (block_support) "block_support" else "random_support")
  .save_and_print_mc(OUT_DIR, stem, num_MC, strat_names, snr_values,
                     fdr_mat, tpr_mat, avg_L_mat, avg_T_mat)
  cat("\n")
}


# ==============================================================================
# Main
# --------------
# Run demo 04
# ==============================================================================
if (run_trex_04) {
  trx_04_lloop_strategies(NUM_MC, rnd_coef = FALSE, block_support = FALSE)
  # trx_04_lloop_strategies(NUM_MC, rnd_coef = FALSE, block_support = TRUE)
}
cat("\ndemo_trex_04_lloop_strategies.R complete.\n")
