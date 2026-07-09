# =================================================================================
# demo_trex_gvs_07.R
# =================================================================================
#
# T-Rex+GVS on the ARMA mixed-structure blocks DGP.
#
# Four blocks (sizes {20, 50, 80, 65}) each following a different time-series
# model:
#   Block 1 = AR(1)      [size 20]  — active
#   Block 2 = MA(3)      [size 50]  — active
#   Block 3 = ARMA(2,1)  [size 80]  — active   (s = 150)
#   Block 4 = AR(1) trap [size 65]  — inactive
# Shuffled into p = 500 columns with white-noise gaps. The AR component is
# governed by ar_coef; ma_coefs = c(0.5, 0.3, 0.1) is fixed. This is the most
# heterogeneous within-block structure in the suite.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
# The swept correlation knob here is ar_coef (not rho).
#
# Parts:
#   Part 1: Single-run demo — MA(3) lag-structure check + block layout + one run.
#   Part 2: SNR sweep       — fixed ar_coef = 0.8.
#   Part 3: ar_coef sweep   — fixed SNR = 2.0.
#   Part 4: 2-D SNR x ar_coef grid.
#
# =================================================================================

library(TRexSelectorNeo)
library(plotly)
library(parallel)

num_cores <- local({
  a <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1L]))
  if (!is.na(a) && a > 0L) a else 6L
})

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    flag <- grep("^--file=", args, value = TRUE)
    if (length(flag) > 0L)
      dirname(normalizePath(sub("^--file=", "", flag[1L])))
    else "."
  }
)

source(file.path(this_dir_, "..", "..", "support_generators.R"))
source(file.path(this_dir_, "..", "..", "simulation_utils.R"))
source(file.path(this_dir_, "..", "trex_gvs_dgps.R"))


# ==============================================================================
# Configuration
# ==============================================================================

CFG <- list(
  n        = 200L,
  p        = 500L,
  ar_coef  = 0.8,          # fixed AR parameter for the SNR sweep
  snr      = 2.0,          # fixed SNR for the ar_coef sweep
  K        = 20L,
  tFDR     = 0.1,
  seed     = 2026L,
  num_MC   = 200L,
  corr_max = 0.98,
  hc_dist  = "single"
)

GVS_METHODS <- list(
  list(label = "EN",     GVS_type = "EN",  en_solver = "TENET"),
  list(label = "EN+AUG", GVS_type = "EN",  en_solver = "TENET_AUG"),
  list(label = "IEN",    GVS_type = "IEN", en_solver = "TENET")
)

run_part_1 <- TRUE
run_part_2 <- TRUE
run_part_3 <- TRUE
run_part_4 <- TRUE


# ==============================================================================
# Local helpers
# ==============================================================================

#' Print a formatted single-run result block.
.print_result <- function(scenario_name, dat, result, tFDR) {
  ts      <- which(dat$beta != 0)
  sel_set <- result$selected_indices
  n_sel   <- length(sel_set)
  n_tp    <- length(intersect(sel_set, ts))
  n_fp    <- n_sel - n_tp

  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf(" Data: n = %d, p = %d, tFDR = %.2f, s = %d\n",
              dat$n, dat$p, tFDR, dat$s))
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f\n", dat$snr, dat$sigma_y))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Calibration:  T_stop = %d,  L = %d\n",
              result$T_stop, result$L))
  cat(sprintf("  Selection:    %d selected  |  TP = %d  FP = %d\n",
              n_sel, n_tp, n_fp))
  cat(sprintf("  Rates:        TPP = %.3f  |  FDP = %.3f  (target tFDR <= %.2f)\n",
              TRexSelectorNeo::compute_tpp(sel_set, ts),
              TRexSelectorNeo::compute_fdp(sel_set, ts), tFDR))
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}

#' Run one MC point on the ARMA-blocks DGP for a given (method, snr, ar_coef).
.arma_point <- function(method, snr, ar_coef) {
  .run_mc_arma_blocks(
    n = CFG$n, p = CFG$p, snr = snr, ar_coef = ar_coef,
    tFDR = CFG$tFDR, K = CFG$K, num_MC = CFG$num_MC, seed = CFG$seed,
    GVS_type = method$GVS_type, corr_max = CFG$corr_max, hc_dist = CFG$hc_dist,
    en_solver = method$en_solver, n_cores = num_cores)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("ARMA-Blocks GVS  |  Part 1: Single-run\n")
  cat(sprintf("n=%d,  p=%d,  s=150  (blocks 20+50+80 active; 65 trap)\n",
              CFG$n, CFG$p))
  cat(sprintf("ar_coef=%.2f,  SNR=%.2f,  tFDR=%.2f\n",
              CFG$ar_coef, CFG$snr, CFG$tFDR))
  cat("Block types: AR(1) | MA(3) | ARMA(2,1) | AR(1) trap\n")
  cat(strrep("=", 70), "\n\n")

  cat("[Part 1] Generating ARMA-blocks data ...\n")
  dat_p1 <- dgp_arma_blocks_snr(n = CFG$n, p = CFG$p, snr = CFG$snr,
                                ar_coef = CFG$ar_coef, seed = CFG$seed)
  cat(sprintf("[Part 1] Active: %d  |  sigma_y = %.4f\n\n",
              dat_p1$s, dat_p1$sigma_y))

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

  # MA(3) lag structure: positive short-lag correlation, ~zero beyond lag 3.
  ma3_idx <- dat_p1$block_indices[[2L]]
  lag1 <- cor(dat_p1$X[, ma3_idx[1L]], dat_p1$X[, ma3_idx[2L]])
  lag4 <- cor(dat_p1$X[, ma3_idx[1L]], dat_p1$X[, ma3_idx[5L]])
  cat(sprintf("[Part 1] MA(3) Block 2 — Lag 1 correlation (expect > 0): %.4f\n", lag1))
  cat(sprintf("[Part 1] MA(3) Block 2 — Lag 4 (expect ~0):              %.4f\n\n", lag4))

  cormat_p1 <- plot_cormat(cor(dat_p1$X))
  if (interactive()) {
    cat("[Part 1] Displaying correlation matrix ...\n")
    print(cormat_p1)
  }

  for (m in GVS_METHODS) {
    cat(sprintf("[Part 1] Running trex+GVS (%s) ...\n\n", m$label))
    res <- .gvs_select(dat_p1$X, dat_p1$y, CFG$tFDR, CFG$K,
                       m$GVS_type, CFG$corr_max, CFG$hc_dist, m$en_solver)
    .print_result(sprintf("Part 1 — ARMA-Blocks GVS  [%s]", m$label),
                  dat_p1, res, CFG$tFDR)
  }
}


# ==============================================================================
# Part 2: MC simulation - SNR sweep
# ==============================================================================
# Fixed ar_coef = 0.8. Compares EN / EN+AUG / IEN.

if (run_part_2) {
  snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

  cat(strrep("=", 70), "\n")
  cat("ARMA-Blocks GVS MC — Part 2: SNR sweep\n")
  cat(sprintf("tFDR=%.2f  ar_coef=%.2f  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$ar_coef, CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] SNR sweep, %d MC per point ...\n", m$label, CFG$num_MC))
    res <- lapply(snr_grid, function(s) {
      cat(sprintf("    SNR = %.2f ...\n", s))
      c(list(SNR = s), .arma_point(m, s, CFG$ar_coef))
    })
    cat(sprintf("\n  Results — %s:\n", m$label))
    .print_table(res, "SNR")
  }
}


# ==============================================================================
# Part 3: MC simulation - ar_coef sweep
# ==============================================================================
# Fixed SNR = 2.0. ar_coef passed directly to the DGP.

if (run_part_3) {
  ar_coef_grid <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

  cat(strrep("=", 70), "\n")
  cat("ARMA-Blocks GVS MC — Part 3: ar_coef sweep\n")
  cat(sprintf("tFDR=%.2f  SNR=%.2f  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$snr, CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] ar_coef sweep, %d MC per point ...\n",
                m$label, CFG$num_MC))
    res <- lapply(ar_coef_grid, function(ar_val) {
      cat(sprintf("    ar_coef = %.2f ...\n", ar_val))
      c(list(ar_coef = ar_val), .arma_point(m, CFG$snr, ar_val))
    })
    cat(sprintf("\n  Results — %s:\n", m$label))
    .print_table(res, "ar_coef", "%-8.2f")
  }
}


# ==============================================================================
# Part 4: 2-D phase-transition sweep - SNR x ar_coef
# ==============================================================================

if (run_part_4) {
  snr_grid     <- c(0.2, 0.5, 1.0, 2.0, 5.0)
  ar_coef_grid <- c(0.30, 0.50, 0.70, 0.80, 0.90, 0.95)
  row_nms      <- paste0("snr=", sprintf("%.2f", snr_grid))
  col_nms      <- paste0("ar=", sprintf("%.2f", ar_coef_grid))

  cat(strrep("=", 70), "\n")
  cat("ARMA-Blocks GVS MC  |  Part 4: SNR x ar_coef grid\n")
  cat(sprintf("n=%d  p=%d  %d MC  tFDR=%.2f  corr_max=%.2f\n",
              CFG$n, CFG$p, CFG$num_MC, CFG$tFDR, CFG$corr_max))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    mat_TPP <- matrix(NA_real_, length(snr_grid), length(ar_coef_grid),
                      dimnames = list(row_nms, col_nms))
    mat_FDP <- mat_TPP
    for (i in seq_along(snr_grid)) {
      for (j in seq_along(ar_coef_grid)) {
        cat(sprintf("  [%s] SNR=%.2f ar_coef=%.2f ...\n",
                    m$label, snr_grid[i], ar_coef_grid[j]))
        r <- .arma_point(m, snr_grid[i], ar_coef_grid[j])
        mat_TPP[i, j] <- r$mean_TPP
        mat_FDP[i, j] <- r$mean_FDP
      }
    }
    .print_matrix(mat_TPP,
                  sprintf("mean_TPP [%s] (rows: SNR, cols: ar_coef)", m$label))
    .print_matrix(mat_FDP,
                  sprintf("mean_FDP [%s] (rows: SNR, cols: ar_coef)", m$label))
  }
}


cat("\nARMA-blocks GVS demo complete.\n")
