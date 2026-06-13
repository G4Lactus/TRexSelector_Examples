# ==============================================================================
# demo_trex_gvs_05.R
# ==============================================================================
#
# T-Rex+GVS demo and MC simulation for the random-blocks DGP.
# 3 active blocks (sizes 20, 50, 80 = 150 active), shuffled into a random order
# and placed in the p-dimensional space with random white-noise gaps.
#
#  Part 1: Single-run demo — block placement diagnostics + selector run.
#
# Monte Carlo simulations (Parts 2–4) are in simulation_trex_gvs_05.R.
#

# ==============================================================================

library(TRexSelector)
library(plotly)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) "."
)

source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))
source(file.path(this_dir_, "dgp_random_blocks.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n        = 200L,
  p        = 500L,
  sd_x     = sqrt(0.01),
  snr      = 2.0,
  K        = 20L,
  tFDR     = 0.1,
  seed     = 2026L,
  corr_max = 0.98,
  hc_dist  = "single"
)

run_part_1 <- TRUE


# ==============================================================================
# Local helpers
# ==============================================================================

.print_result <- function(scenario_name, dat, result, tFDR) {
  sel_set <- which(result$selected_var == 1L)
  n_sel   <- length(sel_set)
  n_tp    <- length(intersect(sel_set, which(dat$beta != 0)))
  n_fp    <- n_sel - n_tp
  tpp_val <- TRexSelector::TPP(result$selected_var, dat$beta)
  fdp_val <- TRexSelector::FDP(result$selected_var, dat$beta)

  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf(" Data: n = %d, p = %d, tFDR = %.2f, s = %d\n",
              dat$n, dat$p, tFDR, dat$s))
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f,  sd_x = %.4f\n",
              dat$snr, dat$sigma_y, dat$sd_x))
  cat(sprintf("       group_order = {%s}  (sizes={%s})\n",
              paste(dat$group_order, collapse = ", "),
              paste(dat$group_sizes[dat$group_order], collapse = ", ")))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Calibration:  T_stop = %d,  dummies = %d\n",
              result$T_stop, result$num_dummies))
  cat(sprintf("  Selection:    %d selected  |  TP = %d  FP = %d\n",
              n_sel, n_tp, n_fp))
  cat(sprintf("  Rates:        TPP = %.3f  |  FDP = %.3f  (target tFDR <= %.2f)\n",
              tpp_val, fdp_val, tFDR))
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  trx_gvs_05_01 <- function() {

    # Header
    # ----------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Random-Blocks GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  blocks={20,50,80} randomly placed\n",
                PARAMS$n, PARAMS$p))
    cat(sprintf("SNR=%.2f,  sd_x=%.4f,  tFDR=%.2f\n",
                PARAMS$snr, PARAMS$sd_x, PARAMS$tFDR))
    cat(strrep("=", 70), "\n\n")
    # ----------------------------------------------------

    cat("[Part 1] Generating random-blocks data ...\n")
    dat_p1 <- dgp_random_blocks_snr(
      n    = PARAMS$n,
      p    = PARAMS$p,
      snr  = PARAMS$snr,
      sd_x = PARAMS$sd_x,
      seed = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active variables: %d  |  sigma_y = %.4f\n",
                dat_p1$s, dat_p1$sigma_y))
    cat(sprintf("[Part 1] Block placement order (original IDs): {%s}\n",
                paste(dat_p1$group_order, collapse = ", ")))
    for (b in seq_len(3)) {
      idx <- dat_p1$group_indices[[b]]
      cat(sprintf("[Part 1]   Block %d (size %2d): cols %3d \u2013 %3d\n",
                  b, dat_p1$group_sizes[b], min(idx), max(idx)))
    }
    cat("\n")

    cat("[Part 1] Plotting correlation matrix ...\n")
    plot_cormat(cor(dat_p1$X))

    cat("[Part 1] Running trex+GVS (EN) ...\n\n")
    res_en <- TRexSelector::trex(
      X             = dat_p1$X,
      y             = dat_p1$y,
      tFDR          = PARAMS$tFDR,
      K             = PARAMS$K,
      method        = "trex+GVS",
      GVS_type      = "EN",
      lambda_2_lars = NULL,
      hc_dist       = PARAMS$hc_dist,
      corr_max      = PARAMS$corr_max,
      verbose       = FALSE,
      seed          = PARAMS$seed
    )
    .print_result("Part 1 \u2014 Random-Blocks GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)

    cat("[Part 1] Running trex+GVS (IEN) ...\n\n")
    res_ien <- TRexSelector::trex(
      X             = dat_p1$X,
      y             = dat_p1$y,
      tFDR          = PARAMS$tFDR,
      K             = PARAMS$K,
      method        = "trex+GVS",
      GVS_type      = "IEN",
      lambda_2_lars = NULL,
      hc_dist       = PARAMS$hc_dist,
      corr_max      = PARAMS$corr_max,
      verbose       = FALSE,
      seed          = PARAMS$seed
    )
    .print_result("Part 1 \u2014 Random-Blocks GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
  }

  trx_gvs_05_01()

}  # end Part 1
# ==============================================================================



cat("\nRandom-blocks GVS demo complete.\n")
