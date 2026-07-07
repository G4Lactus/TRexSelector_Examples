# ==============================================================================
# simulation_trex_gvs_08.R
# ==============================================================================
#
# T-Rex+GVS Monte Carlo simulations for the ARMA mixed-structure blocks DGP.
# 4 blocks: Block1=AR(1) [size 20], Block2=MA(3) [size 50],
#           Block3=ARMA(2,1) [size 80] -- all active (beta=3);
#           Block4=AR(1) trap [size 65] -- inactive.
# Shuffled into p=500 columns with white-noise gaps. Active set: s = 150.
# ar_coef governs the AR component; ma_coefs = c(0.5, 0.3, 0.1) fixed.
#
#  Part 2: MC simulation sweeping SNR. Compares EN vs IEN.
#          Fixed: ar_coef=0.8, n=200, p=500, tFDR=0.1, K=20, 200 MC.
#          SNR grid: {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
#  Part 3: MC simulation sweeping AR coefficient ar_coef.
#          Fixed: SNR=2.0, n=200, p=500, tFDR=0.1, K=20, 200 MC.
#          ar_coef grid: {0.10, 0.20, 0.30, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95}.
#
#  Part 4: 2D phase-transition sweep: SNR x ar_coef.
#          Grid: SNR in {0.2, 0.5, 1.0, 2.0, 5.0},
#                ar_coef in {0.30, 0.50, 0.70, 0.80, 0.90, 0.95}.
#          Outputs TPP/FDP matrices for EN and IEN.
#
# Single-run demo (Part 1) is in demo_trex_gvs_08.R.
#
# ==============================================================================

run_part_2 <- TRUE
run_part_3 <- TRUE
run_part_4 <- TRUE

library(TRexSelectorNeo)
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

source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))
source(file.path(this_dir_, "dgp_arma_blocks.R"))


# ==============================================================================
# Monte Carlo base parameters
# ==============================================================================

MC_BASE <- list(
  n        = 200L,
  p        = 500L,
  ar_coef  = 0.8,
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
  trx_gvs_11_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    .run_snr_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("ARMA-Blocks GVS MC \u2014 SNR sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  ar_coef: %.2f\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$ar_coef))
      cat(sprintf("n=%d  p=%d  SNR swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(snr_grid, function(snr_val) {
        cat(sprintf("\n  [SNR = %.2f]  running %d MC trials ...\n",
                    snr_val, MC_BASE$num_MC))
        r <- .run_mc_arma_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          ar_coef  = MC_BASE$ar_coef,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = GVS_type,
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist
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

  trx_gvs_11_02()

}  # end Part 2
# ==============================================================================


# ==============================================================================
# Part 3: Monte Carlo simulation — ar_coef sweep
# ==============================================================================

if (run_part_3) {
  trx_gvs_11_03 <- function() {

    ar_coef_grid <- c(0.10, 0.20, 0.30, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)

    .run_ar_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("ARMA-Blocks GVS MC \u2014 ar_coef sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  SNR: %.2f\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$snr))
      cat(sprintf("n=%d  p=%d  ar_coef swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(ar_coef_grid, function(ar_val) {
        cat(sprintf("\n  [ar_coef = %.2f]  running %d MC trials ...\n",
                    ar_val, MC_BASE$num_MC))
        r <- .run_mc_arma_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = MC_BASE$snr,
          ar_coef  = ar_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = GVS_type,
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist
        )
        c(list(ar_coef = ar_val), r)
      })
    }

    res_en <- .run_ar_sweep("EN")
    cat("\n  Results \u2014 EN:\n")
    .print_table(res_en, "ar_coef", "%-8.2f")

    res_ien <- .run_ar_sweep("IEN")
    cat("\n  Results \u2014 IEN:\n")
    .print_table(res_ien, "ar_coef", "%-8.2f")
  }

  trx_gvs_11_03()

}  # end Part 3
# ==============================================================================


# ==============================================================================
# Part 4: 2D phase-transition sweep: SNR x ar_coef
# ==============================================================================

if (run_part_4) {
  trx_gvs_11_04 <- function() {

    snr_grid     <- c(0.2, 0.5, 1.0, 2.0, 5.0)
    ar_coef_grid <- c(0.30, 0.50, 0.70, 0.80, 0.90, 0.95)
    n_snr        <- length(snr_grid)
    n_ar         <- length(ar_coef_grid)

    row_nms <- paste0("snr=", sprintf("%.2f", snr_grid))
    col_nms <- paste0("ar=",  sprintf("%.2f", ar_coef_grid))

    mat_TPP_en  <- matrix(NA_real_, nrow = n_snr, ncol = n_ar,
                          dimnames = list(row_nms, col_nms))
    mat_FDP_en  <- mat_TPP_en
    mat_TPP_ien <- mat_TPP_en
    mat_FDP_ien <- mat_TPP_en

    cat(strrep("=", 70), "\n")
    cat("ARMA-Blocks GVS MC  |  Part 4: SNR x ar_coef sweep\n")
    cat(sprintf("n=%d, p=%d, %d MC, tFDR=%.2f, corr_max=%.2f\n",
                MC_BASE$n, MC_BASE$p, MC_BASE$num_MC,
                MC_BASE$tFDR, MC_BASE$corr_max))
    cat(sprintf("  snr_grid    : {%s}\n", paste(snr_grid, collapse = ", ")))
    cat(sprintf("  ar_coef_grid: {%s}\n",
                paste(sprintf("%.2f", ar_coef_grid), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    total_cells <- n_snr * n_ar
    cell_idx    <- 0L

    for (i_snr in seq_len(n_snr)) {
      snr_val <- snr_grid[i_snr]
      for (i_ar in seq_len(n_ar)) {
        ar_val   <- ar_coef_grid[i_ar]
        cell_idx <- cell_idx + 1L
        prefix <- sprintf("  [%d/%d] SNR=%.2f  ar_coef=%.2f",
                          cell_idx, total_cells, snr_val, ar_val)

        cat(sprintf("%s  [EN]  running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_en <- .run_mc_arma_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          ar_coef  = ar_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = "EN",
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist
        )
        mat_TPP_en[i_snr, i_ar] <- r_en$mean_TPP
        mat_FDP_en[i_snr, i_ar] <- r_en$mean_FDP
        cat(sprintf("    %s [EN]  done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_en$mean_TPP, r_en$mean_FDP))

        cat(sprintf("%s  [IEN] running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_ien <- .run_mc_arma_blocks(
          n        = MC_BASE$n,
          p        = MC_BASE$p,
          snr      = snr_val,
          ar_coef  = ar_val,
          tFDR     = MC_BASE$tFDR,
          K        = MC_BASE$K,
          num_MC   = MC_BASE$num_MC,
          seed     = MC_BASE$seed,
          GVS_type = "IEN",
          corr_max = MC_BASE$corr_max,
          hc_dist  = MC_BASE$hc_dist
        )
        mat_TPP_ien[i_snr, i_ar] <- r_ien$mean_TPP
        mat_FDP_ien[i_snr, i_ar] <- r_ien$mean_FDP
        cat(sprintf("    %s [IEN] done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_ien$mean_TPP, r_ien$mean_FDP))
      }
    }

    .print_matrix(mat_TPP_en,  "mean_TPP [EN]  (rows: SNR, cols: ar_coef)")
    .print_matrix(mat_FDP_en,  "mean_FDP [EN]  (rows: SNR, cols: ar_coef)")
    .print_matrix(mat_TPP_ien, "mean_TPP [IEN] (rows: SNR, cols: ar_coef)")
    .print_matrix(mat_FDP_ien, "mean_FDP [IEN] (rows: SNR, cols: ar_coef)")
  }

  trx_gvs_11_04()

}  # end Part 4
# ==============================================================================


cat("\nARMA-blocks GVS MC simulations complete.\n")
