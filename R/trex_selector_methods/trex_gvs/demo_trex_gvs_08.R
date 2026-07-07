# ==============================================================================
# demo_trex_gvs_08.R
# ==============================================================================
#
# T-Rex+GVS demo and MC simulation for the ARMA mixed-structure blocks DGP.
# 4 blocks: Block1=AR(1) [size 20],
#            Block2=MA(3) [size 50],
#           Block3=ARMA(2,1) [size 80] — all active (beta=3);
#           Block4=AR(1) trap [size 65] — inactive.
# Shuffled into p=500 columns with white-noise gaps. Active set: s = 150.
# ar_coef governs the AR component; ma_coefs = c(0.5, 0.3, 0.1) fixed.
#
#  Part 1: Single-run demo — MA(3) lag structure check + selector.
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
source(file.path(this_dir_, "dgp_arma_blocks.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n        = 200L,
  p        = 500L,
  ar_coef  = 0.8,
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
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f,  ar_coef = %.2f\n",
              dat$snr, dat$sigma_y, dat$ar_coef))
  cat(sprintf("       Blocks (sizes 20,50,80 active; 65 trap), shuffled; ARMA\n"))
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
  trx_gvs_11_01 <- function() {

    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("ARMA-Blocks GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  s=150  (blocks 20+50+80 active; 65 trap)\n",
                PARAMS$n, PARAMS$p))
    cat(sprintf("ar_coef=%.2f,  SNR=%.2f,  tFDR=%.2f\n",
                PARAMS$ar_coef, PARAMS$snr, PARAMS$tFDR))
    cat("Block types: AR(1) | MA(3) | ARMA(2,1) | AR(1) trap\n")
    cat(strrep("=", 70), "\n\n")

    cat("[Part 1] Generating ARMA-blocks data ...\n")
    dat_p1 <- dgp_arma_blocks_snr(
      n       = PARAMS$n,
      p       = PARAMS$p,
      snr     = PARAMS$snr,
      ar_coef = PARAMS$ar_coef,
      seed    = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active: %d  |  sigma_y = %.4f\n\n",
                dat_p1$s, dat_p1$sigma_y))

    # Block layout + order
    cat("[Part 1] Block layout (shuffled order):\n")
    block_sizes  <- c(20L, 50L, 80L, 65L)
    block_types  <- c("AR(1)", "MA(3)", "ARMA(2,1)", "AR(1) trap")
    active_label <- c("active", "active", "active", "trap")
    for (b in seq_len(4L)) {
      idx <- dat_p1$block_indices[[b]]
      cat(sprintf("  Block %d [size=%2d, %s, %s]: columns %d - %d\n",
                  b, block_sizes[b], block_types[b], active_label[b],
                  min(idx), max(idx)))
    }
    cat(sprintf("  Block order (ID sequence placed): {%s}\n\n",
                paste(dat_p1$block_order, collapse = ", ")))

    # MA(3) lag structure check
    ma3_idx <- dat_p1$block_indices[[2L]]
    lag1 <- cor(dat_p1$X[, ma3_idx[1L]], dat_p1$X[, ma3_idx[2L]])
    lag4 <- cor(dat_p1$X[, ma3_idx[1L]], dat_p1$X[, ma3_idx[5L]])
    cat(sprintf("[Part 1] MA(3) Block 2 — Lag 1 correlation (expect > 0): %.4f\n", lag1))
    cat(sprintf("[Part 1] MA(3) Block 2 — Lag 4 (expect ~0):              %.4f\n\n", lag4))

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
    .print_result("Part 1 \u2014 ARMA-Blocks GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)

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
    .print_result("Part 1 \u2014 ARMA-Blocks GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
  }

  trx_gvs_11_01()

}  # end Part 1
# ==============================================================================



cat("\nARMA-blocks GVS demo complete.\n")
