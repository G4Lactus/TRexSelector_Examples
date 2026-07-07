# ==============================================================================
# demo_trex_gvs_04.R
# ==============================================================================
#
# T-Rex+GVS demo and MC simulation for the negative-traps DGP.
# Group 1 (indices 1-100):   active, +Z1/-Z1, beta = +3/-3.  s = 100.
# Group 2 (indices 101-200): inactive Trap 1, +Z2/-Z2.
# Noise   (indices 201-300): white noise.
# Group 3 (indices 301-360): inactive Trap 2, +Z3/-Z3.
# Noise   (indices 361-500): white noise.
#
#  Part 1: Single-run demo — trap FP autopsy + selector run.
#
# Monte Carlo simulations (Parts 2–4) are in simulation_trex_gvs_04.R.
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
source(file.path(this_dir_, "dgp_neg_traps.R"))


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
  cat(sprintf("       Group 1 [1-100]:   active (+Z1/-Z1)\n"))
  cat(sprintf("       Group 2 [101-200]: Trap 1 (+Z2/-Z2)\n"))
  cat(sprintf("       Group 3 [301-360]: Trap 2 (+Z3/-Z3)\n"))
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

#' Categorize false positives by source region (trap 1, trap 2, or noise).
.fp_autopsy <- function(selected_var, beta) {
  sel_set    <- which(selected_var == 1L)
  fp_indices <- sel_set[beta[sel_set] == 0]

  if (length(fp_indices) == 0L) {
    cat("  No false positives.\n\n")
    return(invisible(NULL))
  }

  trap1_leaks <- sum(fp_indices >= 101L & fp_indices <= 200L)
  trap2_leaks <- sum(fp_indices >= 301L & fp_indices <= 360L)
  noise_leaks <- length(fp_indices) - trap1_leaks - trap2_leaks

  cat(sprintf("  --- FP Autopsy (%d total FPs) ---\n", length(fp_indices)))
  cat(sprintf("  Leaked from Trap 1 (101-200): %d\n", trap1_leaks))
  cat(sprintf("  Leaked from Trap 2 (301-360): %d\n", trap2_leaks))
  cat(sprintf("  Leaked from white noise:       %d\n\n", noise_leaks))
  invisible(NULL)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  trx_gvs_07_01 <- function() {

    # Header
    # ----------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Neg-Traps GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  s=100  (active group + 2 inactive traps)\n",
                PARAMS$n, PARAMS$p))
    cat(sprintf("SNR=%.2f,  sd_x=%.4f,  tFDR=%.2f\n",
                PARAMS$snr, PARAMS$sd_x, PARAMS$tFDR))
    cat(strrep("=", 70), "\n\n")
    # ----------------------------------------------------

    cat("[Part 1] Generating neg-traps data ...\n")
    dat_p1 <- dgp_neg_traps_snr(
      n    = PARAMS$n,
      p    = PARAMS$p,
      snr  = PARAMS$snr,
      sd_x = PARAMS$sd_x,
      seed = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active: %d  |  Null (incl. traps): %d  |  sigma_y = %.4f\n\n",
                dat_p1$s, dat_p1$p - dat_p1$s, dat_p1$sigma_y))

    # Structural checks
    cat(sprintf("[Part 1] Cor(X[,1],   X[,51])  [active,  expect ~-%.2f]: %.4f\n",
                1 / (1 + PARAMS$sd_x^2), cor(dat_p1$X[, 1],   dat_p1$X[, 51])))
    cat(sprintf("[Part 1] Cor(X[,301], X[,331]) [trap 2, expect ~-%.2f]: %.4f\n\n",
                1 / (1 + PARAMS$sd_x^2), cor(dat_p1$X[, 301], dat_p1$X[, 331])))

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
    .print_result("Part 1 \u2014 Neg-Traps GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)
    .fp_autopsy(res_en$selected_var, dat_p1$beta)

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
    .print_result("Part 1 \u2014 Neg-Traps GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
    .fp_autopsy(res_ien$selected_var, dat_p1$beta)
  }

  trx_gvs_07_01()

}  # end Part 1
# ==============================================================================



cat("\nNeg-traps GVS demo complete.\n")
