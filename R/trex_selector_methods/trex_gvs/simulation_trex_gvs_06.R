# ==============================================================================
# simulation_trex_gvs_06.R
# ==============================================================================
#
# T-Rex+GVS Monte Carlo simulations for the heavy-tailed equi-correlated
# blocks DGP. Same block structure as demo_08 but with Student-t(3)
# distributed latent factors, noise, and response shocks.
# 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
# Shuffled into p=500 columns with heavy-tailed noise gaps. Active set: s = 150.
#
#  Part 2: MC simulation sweeping SNR. Compares EN vs IEN.
#          Fixed: rho=0.75, df=3, n=200, p=500, tFDR=0.1, K=20, 200 MC.
#          SNR grid: {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
#  Part 3: MC simulation sweeping within-block equi-correlation rho.
#          Fixed: SNR=2.0, n=200, p=500, tFDR=0.1, K=20, 200 MC.
#          rho grid: {0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99}.
#
#  Part 4: 2D phase-transition sweep: SNR x rho.
#          Grid: SNR in {0.2, 0.5, 1.0, 2.0, 5.0},
#                rho in {0.30, 0.50, 0.70, 0.90, 0.95, 0.99}.
#          Outputs TPP/FDP matrices for EN and IEN.
#
# Single-run demo (Part 1) is in demo_trex_gvs_06.R.
#
# ==============================================================================

run_part_2 <- TRUE
run_part_3 <- TRUE
run_part_4 <- TRUE

library(TRexSelector)
library(parallel)

num_cores <- local({
  a <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1L]))
  if (!is.na(a) && a > 0L) a else 6L
})

this_dir_ <- local({
  args <- commandArgs(trailingOnly = FALSE)
  flag <- grep("^--file=", args, value = TRUE)
  if (length(flag) > 0L)
    dirname(normalizePath(sub("^--file=", "", flag[1L])))
  else
    "."
})

log_dir <- file.path(this_dir_, "..", "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
con_log <- file(file.path(log_dir, "simulation_trex_gvs_06.log"), open = "wt")
sink(con_log, split = FALSE)
on.exit({ sink(); close(con_log) }, add = TRUE)

source(file.path(this_dir_, "support_generators.R"))
source(file.path(this_dir_, "simulation_utils.R"))
source(file.path(this_dir_, "dgp_t3_equi_blocks.R"))


# ==============================================================================
# Monte Carlo base parameters
# ==============================================================================

MC_BASE <- list(
  n        = 200L,
  p        = 500L,
  rho      = 0.75,
  df       = 3L,
  snr      = 2.0,
  tFDR     = 0.1,
  K        = 20L,
  num_MC   = 200L,
  seed     = 2026L,
  corr_max = 0.98,
  hc_dist  = "single"
)


# ==============================================================================
# Part 2: Monte Carlo simulation — SNR sweep
# ==============================================================================

if (run_part_2) {
  trx_gvs_09_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    .run_snr_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("HT Equi-Blocks GVS MC \u2014 SNR sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  rho: %.2f  df: %d\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$rho, MC_BASE$df))
      cat(sprintf("n=%d  p=%d  SNR swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(snr_grid, function(snr_val) {
        cat(sprintf("\n  [SNR = %.2f]  running %d MC trials ...\n",
                    snr_val, MC_BASE$num_MC))
        r <- .run_mc_t3_equi_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          rho      = MC_BASE$rho,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = GVS_type,
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist,
          df       = MC_BASE$df
        )
        c(list(SNR = snr_val), r)
      })
    }

    res_en <- .run_snr_sweep("EN")
    cat("\n  Results \u2014 EN:\n")
    .print_table(res_en, "SNR")

    res_ien <- .run_snr_sweep("IEN")
    cat("\n  Results \u2014 IEN:\n")
    .print_table(res_ien, "SNR")
  }

  trx_gvs_09_02()

}  # end Part 2
# ==============================================================================


# ==============================================================================
# Part 3: Monte Carlo simulation — rho sweep
# ==============================================================================

if (run_part_3) {
  trx_gvs_09_03 <- function() {

    rho_grid <- c(0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99)

    .run_rho_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("HT Equi-Blocks GVS MC \u2014 rho sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  SNR: %.2f  df: %d\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$snr, MC_BASE$df))
      cat(sprintf("n=%d  p=%d  rho swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(rho_grid, function(rho_val) {
        cat(sprintf("\n  [rho = %.2f]  running %d MC trials ...\n",
                    rho_val, MC_BASE$num_MC))
        r <- .run_mc_t3_equi_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = MC_BASE$snr,
          rho      = rho_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = GVS_type,
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist,
          df       = MC_BASE$df
        )
        c(list(rho = rho_val), r)
      })
    }

    res_en <- .run_rho_sweep("EN")
    cat("\n  Results \u2014 EN:\n")
    .print_table(res_en, "rho", "%-6.2f")

    res_ien <- .run_rho_sweep("IEN")
    cat("\n  Results \u2014 IEN:\n")
    .print_table(res_ien, "rho", "%-6.2f")
  }

  trx_gvs_09_03()

}  # end Part 3
# ==============================================================================


# ==============================================================================
# Part 4: 2D phase-transition sweep: SNR x rho
# ==============================================================================

if (run_part_4) {
  trx_gvs_09_04 <- function() {

    snr_grid <- c(0.2, 0.5, 1.0, 2.0, 5.0)
    rho_grid <- c(0.30, 0.50, 0.70, 0.90, 0.95, 0.99)
    n_snr    <- length(snr_grid)
    n_rho    <- length(rho_grid)

    row_nms <- paste0("snr=", sprintf("%.2f", snr_grid))
    col_nms <- paste0("rho=", sprintf("%.2f", rho_grid))

    mat_TPP_en  <- matrix(NA_real_, nrow = n_snr, ncol = n_rho,
                          dimnames = list(row_nms, col_nms))
    mat_FDP_en  <- mat_TPP_en
    mat_TPP_ien <- mat_TPP_en
    mat_FDP_ien <- mat_TPP_en

    cat(strrep("=", 70), "\n")
    cat("HT Equi-Blocks GVS MC  |  Part 4: SNR x rho sweep\n")
    cat(sprintf("n=%d, p=%d, df=%d, %d MC, tFDR=%.2f, corr_max=%.2f\n",
                MC_BASE$n, MC_BASE$p, MC_BASE$df, MC_BASE$num_MC,
                MC_BASE$tFDR, MC_BASE$corr_max))
    cat(sprintf("  snr_grid : {%s}\n", paste(snr_grid, collapse = ", ")))
    cat(sprintf("  rho_grid : {%s}\n",
                paste(sprintf("%.2f", rho_grid), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    total_cells <- n_snr * n_rho
    cell_idx    <- 0L

    for (i_snr in seq_len(n_snr)) {
      snr_val <- snr_grid[i_snr]
      for (i_rho in seq_len(n_rho)) {
        rho_val  <- rho_grid[i_rho]
        cell_idx <- cell_idx + 1L
        prefix <- sprintf("  [%d/%d] SNR=%.2f  rho=%.2f",
                          cell_idx, total_cells, snr_val, rho_val)

        cat(sprintf("%s  [EN]  running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_en <- .run_mc_t3_equi_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          rho      = rho_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = "EN",
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist,
          df       = MC_BASE$df
        )
        mat_TPP_en[i_snr, i_rho] <- r_en$mean_TPP
        mat_FDP_en[i_snr, i_rho] <- r_en$mean_FDP
        cat(sprintf("    %s [EN]  done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_en$mean_TPP, r_en$mean_FDP))

        cat(sprintf("%s  [IEN] running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_ien <- .run_mc_t3_equi_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          rho      = rho_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = "IEN",
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist,
          df       = MC_BASE$df
        )
        mat_TPP_ien[i_snr, i_rho] <- r_ien$mean_TPP
        mat_FDP_ien[i_snr, i_rho] <- r_ien$mean_FDP
        cat(sprintf("    %s [IEN] done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_ien$mean_TPP, r_ien$mean_FDP))
      }
    }

    .print_matrix(mat_TPP_en,  "mean_TPP [EN]  (rows: SNR, cols: rho)")
    .print_matrix(mat_FDP_en,  "mean_FDP [EN]  (rows: SNR, cols: rho)")
    .print_matrix(mat_TPP_ien, "mean_TPP [IEN] (rows: SNR, cols: rho)")
    .print_matrix(mat_FDP_ien, "mean_FDP [IEN] (rows: SNR, cols: rho)")
  }

  trx_gvs_09_04()

}  # end Part 4
# ==============================================================================


cat("\nHT equi-blocks GVS MC simulations complete.\n")
