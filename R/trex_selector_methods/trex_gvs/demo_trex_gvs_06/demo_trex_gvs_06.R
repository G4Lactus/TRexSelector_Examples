# =================================================================================
# demo_trex_gvs_06.R
# =================================================================================
#
# T-Rex+GVS on the block-structured AR(1) DGP.
#
# Same 4-block geometry as demo 03 (sizes {20, 50, 80, 65}: 3 active, s = 150,
# plus 1 inactive trap, shuffled into p = 500 with white-noise gaps), but the
# within-block correlation is AR(1)/Toeplitz: Cor(X_j, X_k) = rho^|j-k|, so
# correlation decays with column distance instead of being equicorrelated.
# Robustness to non-equicorrelated within-block structure.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
# Within-block AR(1) parameter rho is passed directly to the DGP.
#
# Parts:
#   Part 1: Single-run demo — AR(1) decay check + block layout + one EN/EN+AUG/IEN run.
#   Part 2: SNR sweep       — fixed rho = 0.85.
#   Part 3: rho sweep       — fixed SNR = 2.0.
#   Part 4: 2-D SNR x rho grid.
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
  rho      = 0.85,         # fixed within-block AR(1) parameter for the SNR sweep
  snr      = 2.0,          # fixed SNR for the rho sweep
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

#' Run one MC point on the AR(1)-blocks DGP for a given (method, snr, rho).
.ar1_point <- function(method, snr, rho) {
  .run_mc_ar1_blocks(
    n = CFG$n, p = CFG$p, snr = snr, rho = rho,
    tFDR = CFG$tFDR, K = CFG$K, num_MC = CFG$num_MC, seed = CFG$seed,
    GVS_type = method$GVS_type, corr_max = CFG$corr_max, hc_dist = CFG$hc_dist,
    en_solver = method$en_solver, n_cores = num_cores)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("AR(1)-Blocks GVS  |  Part 1: Single-run\n")
  cat(sprintf("n=%d,  p=%d,  s=150  (blocks 20+50+80 active; 65 trap)\n",
              CFG$n, CFG$p))
  cat(sprintf("rho=%.2f,  SNR=%.2f,  tFDR=%.2f\n",
              CFG$rho, CFG$snr, CFG$tFDR))
  cat(strrep("=", 70), "\n\n")

  cat("[Part 1] Generating AR(1)-blocks data ...\n")
  dat_p1 <- dgp_ar1_blocks_snr(n = CFG$n, p = CFG$p, snr = CFG$snr,
                               rho = CFG$rho, seed = CFG$seed)
  cat(sprintf("[Part 1] Active: %d  |  sigma_y = %.4f\n\n",
              dat_p1$s, dat_p1$sigma_y))

  cat("[Part 1] Block layout (shuffled order):\n")
  block_sizes  <- c(20L, 50L, 80L, 65L)
  active_label <- c("active", "active", "active", "trap")
  for (b in seq_len(4L)) {
    idx <- dat_p1$block_indices[[b]]
    cat(sprintf("  Block %d [size=%2d, %s]: columns %d - %d\n",
                b, block_sizes[b], active_label[b], min(idx), max(idx)))
  }
  cat(sprintf("  Block order (ID sequence placed): {%s}\n\n",
              paste(dat_p1$block_order, collapse = ", ")))

  # AR(1) decay: empirical vs theoretical rho^k within block 1.
  b1_idx   <- dat_p1$block_indices[[1L]]
  rho_hat1 <- cor(dat_p1$X[, b1_idx[1L]], dat_p1$X[, b1_idx[2L]])
  rho_hat2 <- cor(dat_p1$X[, b1_idx[1L]], dat_p1$X[, b1_idx[3L]])
  rho_hat3 <- cor(dat_p1$X[, b1_idx[1L]], dat_p1$X[, b1_idx[4L]])
  cat("[Part 1] AR(1) decay (theoretical rho^k vs empirical):\n")
  cat(sprintf("  Lag 1: %.4f (target rho   = %.4f)\n", rho_hat1, CFG$rho))
  cat(sprintf("  Lag 2: %.4f (target rho^2 = %.4f)\n", rho_hat2, CFG$rho^2))
  cat(sprintf("  Lag 3: %.4f (target rho^3 = %.4f)\n\n", rho_hat3, CFG$rho^3))

  cormat_p1 <- plot_cormat(cor(dat_p1$X))
  if (interactive()) {
    cat("[Part 1] Displaying correlation matrix ...\n")
    print(cormat_p1)
  }

  for (m in GVS_METHODS) {
    cat(sprintf("[Part 1] Running trex+GVS (%s) ...\n\n", m$label))
    res <- .gvs_select(dat_p1$X, dat_p1$y, CFG$tFDR, CFG$K,
                       m$GVS_type, CFG$corr_max, CFG$hc_dist, m$en_solver)
    .print_result(sprintf("Part 1 — AR(1)-Blocks GVS  [%s]", m$label),
                  dat_p1, res, CFG$tFDR)
  }
}


# ==============================================================================
# Part 2: MC simulation - SNR sweep
# ==============================================================================
# Fixed rho = 0.85. Compares EN / EN+AUG / IEN.

if (run_part_2) {
  snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

  cat(strrep("=", 70), "\n")
  cat("AR(1)-Blocks GVS MC — Part 2: SNR sweep\n")
  cat(sprintf("tFDR=%.2f  rho=%.2f  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$rho, CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] SNR sweep, %d MC per point ...\n", m$label, CFG$num_MC))
    res <- lapply(snr_grid, function(s) {
      cat(sprintf("    SNR = %.2f ...\n", s))
      c(list(SNR = s), .ar1_point(m, s, CFG$rho))
    })
    cat(sprintf("\n  Results — %s:\n", m$label))
    .print_table(res, "SNR")
  }
}


# ==============================================================================
# Part 3: MC simulation - rho sweep
# ==============================================================================
# Fixed SNR = 2.0. rho passed directly to the DGP.

if (run_part_3) {
  rho_grid <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                0.90, 0.95, 0.99)

  cat(strrep("=", 70), "\n")
  cat("AR(1)-Blocks GVS MC — Part 3: rho sweep\n")
  cat(sprintf("tFDR=%.2f  SNR=%.2f  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$snr, CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] rho sweep, %d MC per point ...\n", m$label, CFG$num_MC))
    res <- lapply(rho_grid, function(rho_val) {
      cat(sprintf("    rho = %.2f ...\n", rho_val))
      c(list(rho = rho_val), .ar1_point(m, CFG$snr, rho_val))
    })
    cat(sprintf("\n  Results — %s:\n", m$label))
    .print_table(res, "rho", "%-6.2f")
  }
}


# ==============================================================================
# Part 4: 2-D phase-transition sweep - SNR x rho
# ==============================================================================

if (run_part_4) {
  snr_grid <- c(0.2, 0.5, 1.0, 2.0, 5.0)
  rho_grid <- c(0.30, 0.50, 0.70, 0.85, 0.90, 0.99)
  row_nms  <- paste0("snr=", sprintf("%.2f", snr_grid))
  col_nms  <- paste0("rho=", sprintf("%.2f", rho_grid))

  cat(strrep("=", 70), "\n")
  cat("AR(1)-Blocks GVS MC  |  Part 4: SNR x rho grid\n")
  cat(sprintf("n=%d  p=%d  %d MC  tFDR=%.2f  corr_max=%.2f\n",
              CFG$n, CFG$p, CFG$num_MC, CFG$tFDR, CFG$corr_max))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    mat_TPP <- matrix(NA_real_, length(snr_grid), length(rho_grid),
                      dimnames = list(row_nms, col_nms))
    mat_FDP <- mat_TPP
    for (i in seq_along(snr_grid)) {
      for (j in seq_along(rho_grid)) {
        cat(sprintf("  [%s] SNR=%.2f rho=%.2f ...\n",
                    m$label, snr_grid[i], rho_grid[j]))
        r <- .ar1_point(m, snr_grid[i], rho_grid[j])
        mat_TPP[i, j] <- r$mean_TPP
        mat_FDP[i, j] <- r$mean_FDP
      }
    }
    .print_matrix(mat_TPP, sprintf("mean_TPP [%s] (rows: SNR, cols: rho)", m$label))
    .print_matrix(mat_FDP, sprintf("mean_FDP [%s] (rows: SNR, cols: rho)", m$label))
  }
}


cat("\nAR(1)-blocks GVS demo complete.\n")
