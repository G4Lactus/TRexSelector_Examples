# ==============================================================================
# demo_trex_gvs_12.R
# ==============================================================================
#
# T-Rex+GVS demo and MC simulation for the quasi-hapgen LD-block DGP.
# 7 irregular LD blocks (heterogeneous AR(1)-like within-block correlation) +
# two weak long-range cross-block LD pairs. Fixed p = 500, no block shuffling.
# Active set: blocks 1 (size 30), 3 (size 60), 4 (size 40) -- s = 130.
# Beta = 3 on all active columns.
#
# Block layout (fixed genomic-style order):
#   B1: cols   1-30   rho=0.85  ACTIVE
#   B2: cols  36-115  rho=0.45  trap
#   B3: cols 121-180  rho=0.70  ACTIVE
#   B4: cols 186-225  rho=0.90  ACTIVE
#   B5: cols 231-330  rho=0.55  trap
#   B6: cols 336-385  rho=0.80  trap
#   B7: cols 391-500  rho=0.30  trap (background)
# Long-range LD: cols 10-15 <-> cols 340-345 (str 0.30)
#                cols 130-135 <-> cols 240-245 (str 0.25)
# (Both scaled by rho_scale.)
#
# rho_scale in [0,1] globally multiplies all within-block base correlations
# and long-range LD strengths. Sigma is fixed per rho_scale value, reused
# across all MC trials for reproducibility and speed.
#
#  Part 1: Single-run demo -- LD structure diagnostics + selector.
#
# Monte Carlo simulations (Parts 2–4) are in simulation_trex_gvs_12.R.
#
# ==============================================================================

library(TRexSelector)
library(Matrix)
library(plotly)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) "."
)

source(file.path(this_dir_, "dgp_hapgen_snr.R"))
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n         = 200L,
  p         = 500L,
  s         = 130L,
  rho_scale = 1.0,
  snr       = 2.0,
  K         = 20L,
  tFDR      = 0.1,
  seed      = 2026L,
  corr_max  = 0.98,
  hc_dist   = "single"
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
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f,  rho_scale = %.2f\n",
              dat$snr, dat$sigma_y, dat$rho_scale))
  cat(sprintf("       7 irregular LD blocks; active: B1(30)+B3(60)+B4(40)\n"))
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
  trx_gvs_12_01 <- function() {

    # Header
    # ----------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Quasi-hapgen LD-block GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  s=%d  (B1:30 + B3:60 + B4:40  active)\n",
                PARAMS$n, PARAMS$p, PARAMS$s))
    cat(sprintf("rho_scale=%.2f,  SNR=%.2f,  tFDR=%.2f\n",
                PARAMS$rho_scale, PARAMS$snr, PARAMS$tFDR))
    cat(strrep("=", 70), "\n\n")
    # ----------------------------------------------------

    cat("[Part 1] Generating quasi-hapgen LD-block data ...\n")
    dat_p1 <- dgp_hapgen_snr(
      n         = PARAMS$n,
      snr       = PARAMS$snr,
      rho_scale = PARAMS$rho_scale,
      seed      = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active: %d  |  sigma_y = %.4f\n\n",
                dat_p1$s, dat_p1$sigma_y))

    block_labels <- c("B1 ACTIVE", "B2 trap  ", "B3 ACTIVE", "B4 ACTIVE",
                      "B5 trap  ", "B6 trap  ", "B7 trap  ")
    cat("[Part 1] Fixed block layout:\n")
    for (b in seq_len(7L)) {
      idx <- dat_p1$block_indices[[b]]
      cat(sprintf("  Block %d [%-9s size=%3d  rho_base=%.2f  rho_eff=%.4f]: cols %3d-%3d\n",
                  b, block_labels[b], dat_p1$block_sizes[b],
                  dat_p1$block_rhos_base[b], dat_p1$block_rhos_eff[b],
                  min(idx), max(idx)))
    }
    cat(sprintf("\n  Long-range LD: cols 10-15 <-> cols 340-345 (base strength 0.30)\n"))
    cat(sprintf("                 cols 130-135 <-> cols 240-245 (base strength 0.25)\n\n"))

    b1 <- dat_p1$block_indices[[1L]]
    rho_eff_b1 <- dat_p1$block_rhos_eff[1L]
    cat(sprintf("[Part 1] B1 within-block: Cor(b1[1], b1[2]) = %.4f  (expect rho_eff=%.4f)\n",
                cor(dat_p1$X[, b1[1L]], dat_p1$X[, b1[2L]]), rho_eff_b1))
    cat(sprintf("[Part 1] Long-range LD:   Cor(col 10, col 340) = %.4f  (base strength=%.2f, scaled=%.4f)\n\n",
                cor(dat_p1$X[, 10L], dat_p1$X[, 340L]),
                0.30, min(0.30 * PARAMS$rho_scale, 0.999)))

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
    .print_result("Part 1 - Quasi-hapgen LD-block GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)

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
    .print_result("Part 1 - Quasi-hapgen LD-block GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
  }

  trx_gvs_12_01()

}  # end Part 1
# ==============================================================================



cat("\nQuasi-hapgen LD-block GVS demo complete.\n")
