# ==============================================================================
# demo_trex_scr_01_screen.R
# ==============================================================================
#
# Screen-TRex Selector, part 01: basic in-memory variable screening.
#
# Screen-TRex extends the classical T-Rex Selector with a fast variable-
# *screening* step for ultra-high-dimensional data. It thresholds the same
# dummy-based voting proportion Phi_j that T-Rex accumulates internally:
#   * Ordinary   rule: select {j : Phi_j > 0.5}  (a majority vote).
#   * Bootstrap-CI rule: build a bootstrap confidence band around the estimated
#     FDR curve (R_boot replicates) and pick the threshold from that band.
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/
#     demo_trex_scr_01_screen_trex.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part 1: single-run Ordinary Screen-TRex (high-dimensional, fixed support).
#  Part 2: single-run Bootstrap-CI Screen-TRex (same data).
#  Part 3: Monte Carlo SNR sweep comparing Ordinary vs. Bootstrap-CI screening
#          (FDR, TPR, and the procedure's self-reported Estimated FDR).
#
# The C++ demo uses num_MC = 500; here we use a smaller count (noted below) so
# the demo finishes in ~1 minute. Results are statistically comparable, not
# bit-identical (different RNG streams).
# ==============================================================================

library(TRexSelectorNeo)

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
source(file.path(this_dir_, "trex_scr_sim_utils.R"))

OUT_DIR <- file.path(this_dir_, "simulation_results")

# Run flags
run_part_1 <- TRUE
run_part_2 <- TRUE
run_part_3 <- TRUE

# Fixed true support (C++ 0-based {27,149,398,4,42} -> R 1-based)
TRUE_SUPPORT <- c(28L, 150L, 399L, 5L, 43L)

NUM_MC <- 60L   # C++ uses 500; reduced for a ~1-minute runtime.


# ==============================================================================
# Parts 1 & 2: single-run Ordinary and Bootstrap-CI screening
# ==============================================================================

run_single <- function(bootstrap) {
  n <- 300L; p <- 1000L; snr <- 5.0
  label <- if (bootstrap) "Bootstrap-CI (Phi via bootstrap band)"
           else            "Ordinary (Phi > 0.5)"

  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat(sprintf("Screen-TRex single run  |  %s\n", label))
  cat(sprintf("High-dimensional (p > n): n = %d, p = %d, s = %d, SNR = %.1f\n\n",
              n, p, length(TRUE_SUPPORT), snr))

  cat("Generating synthetic data ...\n")
  dat <- dgp_iid_snr(n, p, TRUE_SUPPORT, coefs = rep(1, length(TRUE_SUPPORT)),
                     snr = snr, seed = 58L)

  cat("Running Screen-TRex selection ...\n")
  sel <- TRexScreeningSelector$new(
    dat$X, dat$y, seed = 42L, verbose = FALSE,
    screen_control = trex_screen_control(trex_method = "TREX",
                                         use_bootstrap_CI = bootstrap),
    control        = make_scr_trex_ctrl()
  )
  sel$select()
  .print_scr_result(sel, dat$true_support)
}


# ==============================================================================
# Part 3: Monte Carlo SNR sweep — Ordinary vs. Bootstrap-CI
# ==============================================================================

run_mc_sweep <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 3: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI)\n")
  cat(sprintf("High-dimensional: n = 300, p = 1000, s = 10, num_MC = %d\n\n", NUM_MC))

  n <- 300L; p <- 1000L; support_size <- 10L
  snr_values <- c(0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
  methods <- default_scr_methods()
  method_names <- vapply(methods, function(m) m$name, character(1))

  init <- function() setNames(
    lapply(method_names, function(x) numeric(length(snr_values))), method_names)
  fdr_map <- init(); tpr_map <- init(); est_fdr_map <- init()

  trex_ctrl <- make_scr_trex_ctrl()

  for (m in methods) {
    cat(strrep("=", 70L), "\n", sep = "")
    cat(sprintf("Method: %s\n", m$name))
    cat(strrep("=", 70L), "\n", sep = "")
    for (si in seq_along(snr_values)) {
      snr <- snr_values[si]
      make_data <- function(seed) {
        set.seed(seed + 500000L)                 # isolated support RNG
        sup <- sample.int(p, support_size)
        dgp_iid_snr(n, p, sup, coefs = rep(1, support_size), snr = snr, seed = seed)
      }
      base_seed <- 24L + (si - 1L) * 1000L
      res <- run_mc_screen(
        NUM_MC, sprintf("%s  SNR=%.2f", m$name, snr), make_data, trex_ctrl,
        screen_args = list(trex_method = m$trex_method, use_bootstrap_CI = m$bootstrap),
        base_seed = base_seed)
      fdr_map[[m$name]][si]     <- res$avg_fdr
      tpr_map[[m$name]][si]     <- res$avg_tpr
      est_fdr_map[[m$name]][si] <- res$avg_est_fdr
    }
  }

  save_and_print_scr_mc(NUM_MC, OUT_DIR, "d01_screen_trex_mc_n300_p1000",
                        snr_values, method_names, fdr_map, tpr_map, est_fdr_map)
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_1) run_single(bootstrap = FALSE)
if (run_part_2) run_single(bootstrap = TRUE)
if (run_part_3) run_mc_sweep()
cat("\ndemo_trex_scr_01_screen.R complete.\n")
