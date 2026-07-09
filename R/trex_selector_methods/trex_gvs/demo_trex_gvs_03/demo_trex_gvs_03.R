# =================================================================================
# demo_trex_gvs_03.R
# =================================================================================
#
# T-Rex+GVS on the mixed-blocks DGP.
#
# Four unequal blocks of sizes {20, 50, 80, 65}: three active (beta = 3, s = 150)
# plus one inactive equicorrelated trap block, placed in a random order separated
# by white-noise gaps. A more realistic grouped design than the contiguous
# Hastie blocks of demo 01.
#
# Methods compared: EN (TENET), EN+AUG (TENET_AUG), IEN.
#
# Following cpp demo 03, two presets are swept (each = SNR + rho + 2-D grid):
#   * Mixed-Blocks        — default within-block noise sd_x = sqrt(0.01) (rho ~ 0.99).
#   * Mixed-Blocks-Rho075 — fixed rho = 0.75 for the SNR sweep; a wider rho grid.
#
# Parts:
#   Part 1: Single-run demo — block-layout diagnostics + one EN/EN+AUG/IEN run.
#   Part 2: Mixed-Blocks preset        — SNR sweep + rho sweep + 2-D grid.
#   Part 3: Mixed-Blocks-Rho075 preset — SNR sweep + rho sweep + 2-D grid.
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
  sd_x     = sqrt(0.01),   # default within-block noise sd (rho ~ 0.99)
  snr      = 2.0,          # fixed SNR for the rho sweeps
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

# Two presets swept in Parts 2 and 3 (mirrors cpp demo 03).
PRESETS <- list(
  default_sd_x = list(
    tag          = "Mixed-Blocks",
    snr_fixed_sdx = sqrt(0.01),                     # rho ~ 0.99
    snr_note     = "sd_x=sqrt(0.01)",
    snr_grid     = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
    rho_grid     = c(0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99),
    snr_grid_2d  = c(0.2, 0.5, 1.0, 2.0, 5.0),
    rho_grid_2d  = c(0.30, 0.50, 0.70, 0.90, 0.95, 0.99)
  ),
  fixed_rho_075 = list(
    tag          = "Mixed-Blocks-Rho075",
    snr_fixed_sdx = sqrt((1 - 0.75) / 0.75),        # rho = 0.75
    snr_note     = "rho=0.75",
    snr_grid     = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
    rho_grid     = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
                     0.70, 0.80, 0.90, 0.95, 0.99),
    snr_grid_2d  = c(0.2, 0.5, 1.0, 2.0, 5.0),
    rho_grid_2d  = c(0.30, 0.50, 0.70, 0.85, 0.90, 0.99)
  )
)

run_part_1 <- TRUE
run_part_2 <- TRUE
run_part_3 <- TRUE


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

#' Run one MC point on the mixed-blocks DGP for a given (method, snr, sd_x).
.mixed_point <- function(method, snr, sd_x) {
  .run_mc_mixed_blocks(
    n = CFG$n, p = CFG$p, snr = snr, sd_x = sd_x,
    tFDR = CFG$tFDR, K = CFG$K, num_MC = CFG$num_MC, seed = CFG$seed,
    GVS_type = method$GVS_type, corr_max = CFG$corr_max, hc_dist = CFG$hc_dist,
    en_solver = method$en_solver, n_cores = num_cores)
}

#' Run the SNR + rho + 2-D benchmark for one preset.
.run_mixed_benchmark <- function(preset) {
  cat(strrep("=", 70), "\n")
  cat(sprintf("Mixed-Blocks GVS MC  |  Preset: %s\n", preset$tag))
  cat(sprintf("n=%d  p=%d  %d MC  tFDR=%.2f  corr_max=%.2f\n",
              CFG$n, CFG$p, CFG$num_MC, CFG$tFDR, CFG$corr_max))
  cat(strrep("=", 70), "\n\n")

  # --- SNR sweep (fixed sd_x) ---
  cat(sprintf("  [%s] SNR sweep  (fixed %s)\n", preset$tag, preset$snr_note))
  for (m in GVS_METHODS) {
    res <- lapply(preset$snr_grid, function(s) {
      cat(sprintf("    [%s] SNR = %.2f ...\n", m$label, s))
      c(list(SNR = s), .mixed_point(m, s, preset$snr_fixed_sdx))
    })
    cat(sprintf("\n  SNR sweep — %s:\n", m$label))
    .print_table(res, "SNR")
  }

  # --- rho sweep (fixed SNR) ---
  cat(sprintf("  [%s] rho sweep  (fixed SNR = %.2f)\n", preset$tag, CFG$snr))
  for (m in GVS_METHODS) {
    res <- lapply(preset$rho_grid, function(rho_val) {
      sd_x_val <- sqrt((1 - rho_val) / rho_val)
      cat(sprintf("    [%s] rho = %.2f  (sd_x = %.4f) ...\n",
                  m$label, rho_val, sd_x_val))
      c(list(rho = rho_val), .mixed_point(m, CFG$snr, sd_x_val))
    })
    cat(sprintf("\n  rho sweep — %s:\n", m$label))
    .print_table(res, "rho", "%-6.2f")
  }

  # --- 2-D SNR x rho grid ---
  cat(sprintf("  [%s] 2-D SNR x rho grid\n", preset$tag))
  row_nms <- paste0("snr=", sprintf("%.2f", preset$snr_grid_2d))
  col_nms <- paste0("rho=", sprintf("%.2f", preset$rho_grid_2d))
  for (m in GVS_METHODS) {
    mat_TPP <- matrix(NA_real_, length(preset$snr_grid_2d),
                      length(preset$rho_grid_2d),
                      dimnames = list(row_nms, col_nms))
    mat_FDP <- mat_TPP
    for (i in seq_along(preset$snr_grid_2d)) {
      for (j in seq_along(preset$rho_grid_2d)) {
        sd_x_val <- sqrt((1 - preset$rho_grid_2d[j]) / preset$rho_grid_2d[j])
        cat(sprintf("    [%s] SNR=%.2f rho=%.2f ...\n",
                    m$label, preset$snr_grid_2d[i], preset$rho_grid_2d[j]))
        r <- .mixed_point(m, preset$snr_grid_2d[i], sd_x_val)
        mat_TPP[i, j] <- r$mean_TPP
        mat_FDP[i, j] <- r$mean_FDP
      }
    }
    .print_matrix(mat_TPP, sprintf("mean_TPP [%s] (rows: SNR, cols: rho)", m$label))
    .print_matrix(mat_FDP, sprintf("mean_FDP [%s] (rows: SNR, cols: rho)", m$label))
  }
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Mixed-Blocks GVS  |  Part 1: Single-run\n")
  cat(sprintf("n=%d,  p=%d,  active blocks={20,50,80}, inactive={65}\n",
              CFG$n, CFG$p))
  cat(sprintf("SNR=%.2f,  sd_x=%.4f,  tFDR=%.2f\n",
              CFG$snr, CFG$sd_x, CFG$tFDR))
  cat(strrep("=", 70), "\n\n")

  cat("[Part 1] Generating mixed-blocks data ...\n")
  dat_p1 <- dgp_mixed_blocks_snr(n = CFG$n, p = CFG$p, snr = CFG$snr,
                                 sd_x = CFG$sd_x, seed = CFG$seed)
  cat(sprintf("[Part 1] Active variables: %d  |  sigma_y = %.4f\n",
              dat_p1$s, dat_p1$sigma_y))
  cat(sprintf("[Part 1] Block order (IDs): {%s}\n",
              paste(dat_p1$block_order, collapse = ", ")))
  for (b in seq_len(4)) {
    idx <- dat_p1$block_indices[[b]]
    lbl <- if (dat_p1$is_active[b]) "ACTIVE  " else "inactive"
    cat(sprintf("[Part 1]   Block %d (size %2d, %s): cols %3d – %3d\n",
                b, dat_p1$block_sizes[b], lbl, min(idx), max(idx)))
  }
  cat("\n")

  cormat_p1 <- plot_cormat(cor(dat_p1$X))
  if (interactive()) {
    cat("[Part 1] Displaying correlation matrix ...\n")
    print(cormat_p1)
  }

  for (m in GVS_METHODS) {
    cat(sprintf("[Part 1] Running trex+GVS (%s) ...\n\n", m$label))
    res <- .gvs_select(dat_p1$X, dat_p1$y, CFG$tFDR, CFG$K,
                       m$GVS_type, CFG$corr_max, CFG$hc_dist, m$en_solver)
    .print_result(sprintf("Part 1 — Mixed-Blocks GVS  [%s]", m$label),
                  dat_p1, res, CFG$tFDR)
  }
}


# ==============================================================================
# Part 2: Mixed-Blocks preset (default sd_x)
# ==============================================================================

if (run_part_2) .run_mixed_benchmark(PRESETS$default_sd_x)


# ==============================================================================
# Part 3: Mixed-Blocks-Rho075 preset (fixed rho = 0.75)
# ==============================================================================

if (run_part_3) .run_mixed_benchmark(PRESETS$fixed_rho_075)


cat("\nMixed-blocks GVS demo complete.\n")
