# ==============================================================================
# demo_trex_01_single_run.R
# ==============================================================================
#
# Classical T-Rex Selector demo, part 01 is a single-run demonstration of the
# with the TLARS solver.
# A data-generating process functionality `dgp_gauss_snr` supports the demo
# and allows for modifications.
#
# The file mirrors:
# cpp/trex_selector_methods/trex/demo_trex_01_single_run.cpp
#
# Demo content:
# ----------------------------------------------------------------
#
#  - Low- and high-dimensional settings.
#  - Fixed true support; TLARS solver.
#  - Illustrates the reporting of TRex internal statistics
#    phi_prime, phi_mat, fdp_hat_mat, r_mat, voting_grid.
#
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
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "trex_sim_utils.R"))

# ==============================================================================
# Global Simulation Parameters
# ==============================================================================

# Fixed true support (0-based in C++, converted to 1-based)
# C++ (0-based): {27, 149, 43, 128, 42, 4}
# R (1-based): c(28, 150, 44, 129, 43, 5)
TRUE_SUPPORT_P1 <- c(28L, 150L, 44L, 129L, 43L, 5L)

run_part_01 <- TRUE


# ==============================================================================
# Part 1: Single-run basic demo
# ==============================================================================

trx_01_part_01 <- function() {

  # Parameters matching C++ demo_TRexSelector()
  tFDR <- 0.1
  snr  <- 1.0

  ctrl <- trex_control(
    solver                   = "TLARS",
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
    seed <- 4231

    cat("\n", strrep("=", 70L), "\n", sep = "")
    cat("Part 1: Single-run basic T-Rex demo  |  ", dim_label, "\n", sep = "")
    cat(sprintf("n = %d,  p = %d,  tFDR = %.2f,  SNR = %.1f\n\n", n, p, tFDR, snr))

    cat("Generating synthetic data ...\n")
    dat <- dgp_gauss_snr(n, p, TRUE_SUPPORT_P1, amplitude = 1, snr = snr, seed = seed)

    cat("Running T-Rex Selector (TLARS) ...\n\n")
    selector <- TRexSelector$new(
      X       = dat$X,
      y       = dat$y,
      tFDR    = tFDR,
      seed    = seed, # using a != -1 seed for reproducibility in the demo
      verbose = TRUE,
      control = ctrl
    )
    selector$select()

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
# Main
# ------------
# Run the demo
# ==============================================================================

if (run_part_01) trx_01_part_01()
cat("\ndemo_trex_01_single_run.R complete.\n")
