# ==============================================================================
# demo_trex_gvs_08.R
# ==============================================================================
#
# T-Rex+GVS demo and MC simulation for the equi-correlated blocks DGP.
# 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
# Shuffled into p=500 columns with white-noise gaps. Active set: s = 150.
# Within-block equi-correlation: Cor(X_j, X_k) = rho (factor-model construction).
#
#  Part 1: Single-run demo — block layout + correlation check + selector.
#
# Monte Carlo simulations (Parts 2–4) are in simulation_trex_gvs_08.R.
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
source(file.path(this_dir_, "dgp_equi_blocks.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n        = 200L,
  p        = 500L,
  rho      = 0.75,
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
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f,  rho = %.2f\n",
              dat$snr, dat$sigma_y, dat$rho))
  cat(sprintf("       Blocks (sizes 20,50,80 active; 65 trap), shuffled\n"))
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
  trx_gvs_08_01 <- function() {

    # Header
    # ----------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Equi-Blocks GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  s=150  (blocks 20+50+80 active; 65 trap)\n",
                PARAMS$n, PARAMS$p))
    cat(sprintf("rho=%.2f,  SNR=%.2f,  tFDR=%.2f\n",
                PARAMS$rho, PARAMS$snr, PARAMS$tFDR))
    cat(strrep("=", 70), "\n\n")
    # ----------------------------------------------------

    cat("[Part 1] Generating equi-correlated blocks data ...\n")
    dat_p1 <- dgp_equi_blocks_snr(
      n    = PARAMS$n,
      p    = PARAMS$p,
      snr  = PARAMS$snr,
      rho  = PARAMS$rho,
      seed = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active: %d  |  sigma_y = %.4f\n\n",
                dat_p1$s, dat_p1$sigma_y))

    # Block layout
    cat("[Part 1] Block layout (shuffled order):\n")
    block_sizes  <- c(20L, 50L, 80L, 65L)
    active_label <- c("active", "active", "active", "trap")
    for (b in seq_len(4L)) {
      idx <- dat_p1$block_indices[[b]]
      cat(sprintf("  Block %d [size=%2d, %s]: columns %d - %d\n",
                  b, block_sizes[b], active_label[b],
                  min(idx), max(idx)))
    }
    cat(sprintf("  Block order (ID sequence placed): {%s}\n\n",
                paste(dat_p1$block_order, collapse = ", ")))

    # Correlation check
    b1_idx <- dat_p1$block_indices[[1L]]
    cat(sprintf("[Part 1] Cor(X[,b1[1]], X[,b1[2]]) [expect rho=%.2f]: %.4f\n\n",
                PARAMS$rho, cor(dat_p1$X[, b1_idx[1L]], dat_p1$X[, b1_idx[2L]])))

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
    .print_result("Part 1 \u2014 Equi-Blocks GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)

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
    .print_result("Part 1 \u2014 Equi-Blocks GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
  }

  trx_gvs_08_01()

}  # end Part 1
# ==============================================================================



cat("\nEqui-blocks GVS demo complete.\n")
