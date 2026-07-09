# =================================================================================
# demo_trex_gvs_04.R
# =================================================================================
#
# T-Rex+GVS on the negative-traps DGP.
#
#   Group 1 (cols 1-100):   ACTIVE, +Z1 / -Z1, beta = +3 / -3.  s = 100.
#   Group 2 (cols 101-200): inactive Trap 1, +Z2 / -Z2.
#   Noise   (cols 201-300): white noise.
#   Group 3 (cols 301-360): inactive Trap 2, +Z3 / -Z3.
#   Noise   (cols 361-500): white noise.
#
# The active group is negatively correlated within itself (sign-flipped halves)
# and two spatially separated sign-flipped trap groups tempt false positives.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
#
# Parts:
#   Part 1: Single-run demo — negative-correlation checks, one EN/EN+AUG/IEN run,
#                             and a false-positive autopsy by trap.
#   Part 2: SNR sweep       — fixed sd_x = sqrt(0.01) (rho ~ 0.99).
#   Part 3: rho sweep       — fixed SNR = 2.0 (sd_x derived from rho).
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
  sd_x     = sqrt(0.01),   # rho ~ 0.99
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
  cat(sprintf("       SNR = %.2f,  sigma_y = %.4f,  sd_x = %.4f\n",
              dat$snr, dat$sigma_y, dat$sd_x))
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

#' Classify false positives by their origin (trap 1, trap 2, or white noise).
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
  cat(sprintf("  Leaked from white noise:      %d\n\n", noise_leaks))
  invisible(NULL)
}

#' Run one MC point on the neg-traps DGP for a given (method, snr, sd_x).
.neg_traps_point <- function(method, snr, sd_x) {
  .run_mc_neg_traps(
    n = CFG$n, p = CFG$p, snr = snr, sd_x = sd_x,
    tFDR = CFG$tFDR, K = CFG$K, num_MC = CFG$num_MC, seed = CFG$seed,
    GVS_type = method$GVS_type, corr_max = CFG$corr_max, hc_dist = CFG$hc_dist,
    en_solver = method$en_solver, n_cores = num_cores)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Neg-Traps GVS  |  Part 1: Single-run\n")
  cat(sprintf("n=%d,  p=%d,  s=100  (active group + 2 inactive traps)\n",
              CFG$n, CFG$p))
  cat(sprintf("SNR=%.2f,  sd_x=%.4f,  tFDR=%.2f\n",
              CFG$snr, CFG$sd_x, CFG$tFDR))
  cat(strrep("=", 70), "\n\n")

  cat("[Part 1] Generating neg-traps data ...\n")
  dat_p1 <- dgp_neg_traps_snr(n = CFG$n, p = CFG$p, snr = CFG$snr,
                              sd_x = CFG$sd_x, seed = CFG$seed)
  cat(sprintf("[Part 1] Active: %d  |  Null (incl. traps): %d  |  sigma_y = %.4f\n\n",
              dat_p1$s, dat_p1$p - dat_p1$s, dat_p1$sigma_y))

  # Structural checks: sign-flipped halves induce negative within-group correlation.
  cat(sprintf("[Part 1] Cor(X[,1],   X[,51])  [active, expect ~-%.2f]: %.4f\n",
              1 / (1 + CFG$sd_x^2), cor(dat_p1$X[, 1],   dat_p1$X[, 51])))
  cat(sprintf("[Part 1] Cor(X[,301], X[,331]) [trap 2, expect ~-%.2f]: %.4f\n\n",
              1 / (1 + CFG$sd_x^2), cor(dat_p1$X[, 301], dat_p1$X[, 331])))

  cormat_p1 <- plot_cormat(cor(dat_p1$X))
  if (interactive()) {
    cat("[Part 1] Displaying correlation matrix ...\n")
    print(cormat_p1)
  }

  for (m in GVS_METHODS) {
    cat(sprintf("[Part 1] Running trex+GVS (%s) ...\n\n", m$label))
    res <- .gvs_select(dat_p1$X, dat_p1$y, CFG$tFDR, CFG$K,
                       m$GVS_type, CFG$corr_max, CFG$hc_dist, m$en_solver)
    .print_result(sprintf("Part 1 — Neg-Traps GVS  [%s]", m$label),
                  dat_p1, res, CFG$tFDR)
    .fp_autopsy(res$selected_var, dat_p1$beta)
  }
}


# ==============================================================================
# Part 2: MC simulation - SNR sweep
# ==============================================================================
# Fixed sd_x = sqrt(0.01) (rho ~ 0.99). Compares EN / EN+AUG / IEN.

if (run_part_2) {
  snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

  cat(strrep("=", 70), "\n")
  cat("Neg-Traps GVS MC — Part 2: SNR sweep\n")
  cat(sprintf("tFDR=%.2f  sd_x=%.4f (rho~%.2f)  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$sd_x, 1 / (1 + CFG$sd_x^2),
              CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] SNR sweep, %d MC per point ...\n", m$label, CFG$num_MC))
    res <- lapply(snr_grid, function(s) {
      cat(sprintf("    SNR = %.2f ...\n", s))
      c(list(SNR = s), .neg_traps_point(m, s, CFG$sd_x))
    })
    cat(sprintf("\n  Results — %s:\n", m$label))
    .print_table(res, "SNR")
  }
}


# ==============================================================================
# Part 3: MC simulation - rho sweep
# ==============================================================================
# Fixed SNR = 2.0. sd_x derived as sqrt((1 - rho) / rho).

if (run_part_3) {
  rho_grid <- c(0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99)

  cat(strrep("=", 70), "\n")
  cat("Neg-Traps GVS MC — Part 3: rho sweep\n")
  cat(sprintf("tFDR=%.2f  SNR=%.2f  n=%d  p=%d  %d MC\n",
              CFG$tFDR, CFG$snr, CFG$n, CFG$p, CFG$num_MC))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    cat(sprintf("\n  [%s] rho sweep, %d MC per point ...\n", m$label, CFG$num_MC))
    res <- lapply(rho_grid, function(rho_val) {
      sd_x_val <- sqrt((1 - rho_val) / rho_val)
      cat(sprintf("    rho = %.2f  (sd_x = %.4f) ...\n", rho_val, sd_x_val))
      c(list(rho = rho_val), .neg_traps_point(m, CFG$snr, sd_x_val))
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
  rho_grid <- c(0.30, 0.50, 0.70, 0.90, 0.95, 0.99)
  row_nms  <- paste0("snr=", sprintf("%.2f", snr_grid))
  col_nms  <- paste0("rho=", sprintf("%.2f", rho_grid))

  cat(strrep("=", 70), "\n")
  cat("Neg-Traps GVS MC  |  Part 4: SNR x rho grid\n")
  cat(sprintf("n=%d  p=%d  %d MC  tFDR=%.2f  corr_max=%.2f\n",
              CFG$n, CFG$p, CFG$num_MC, CFG$tFDR, CFG$corr_max))
  cat(strrep("=", 70), "\n\n")

  for (m in GVS_METHODS) {
    mat_TPP <- matrix(NA_real_, length(snr_grid), length(rho_grid),
                      dimnames = list(row_nms, col_nms))
    mat_FDP <- mat_TPP
    for (i in seq_along(snr_grid)) {
      for (j in seq_along(rho_grid)) {
        sd_x_val <- sqrt((1 - rho_grid[j]) / rho_grid[j])
        cat(sprintf("  [%s] SNR=%.2f rho=%.2f ...\n",
                    m$label, snr_grid[i], rho_grid[j]))
        r <- .neg_traps_point(m, snr_grid[i], sd_x_val)
        mat_TPP[i, j] <- r$mean_TPP
        mat_FDP[i, j] <- r$mean_FDP
      }
    }
    .print_matrix(mat_TPP, sprintf("mean_TPP [%s] (rows: SNR, cols: rho)", m$label))
    .print_matrix(mat_FDP, sprintf("mean_FDP [%s] (rows: SNR, cols: rho)", m$label))
  }
}


cat("\nNeg-traps GVS demo complete.\n")
