# ==============================================================================
# simulation_trex_gvs_12.R
# ==============================================================================
#
# T-Rex+GVS Monte Carlo simulations for the quasi-hapgen LD-block DGP.
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
#  Part 2: MC simulation sweeping SNR. Compares EN vs IEN.
#          Fixed: rho_scale=1.0, n=200, tFDR=0.1, K=20, 200 MC.
#          SNR grid: {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
#  Part 3: MC simulation sweeping rho_scale.
#          Fixed: SNR=2.0, n=200, tFDR=0.1, K=20, 200 MC.
#          rho_scale grid: {0.15, 0.30, 0.50, 0.70, 0.85, 1.00}.
#
#  Part 4: 2D phase-transition sweep: SNR x rho_scale.
#          Grid: SNR in {0.2, 0.5, 1.0, 2.0, 5.0},
#                rho_scale in {0.30, 0.50, 0.70, 0.85, 1.00}.
#          Outputs TPP/FDP matrices for EN and IEN.
#
# Single-run demo (Part 1) is in demo_trex_gvs_12.R.
#
# ==============================================================================

run_part_2 <- TRUE
run_part_3 <- TRUE
run_part_4 <- TRUE

library(TRexSelector)
library(Matrix)
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
con_log <- file(file.path(log_dir, "simulation_trex_gvs_12.log"), open = "wt")
sink(con_log, split = FALSE)
on.exit({ sink(); close(con_log) }, add = TRUE)

source(file.path(this_dir_, "dgp_hapgen_snr.R"))
source(file.path(this_dir_, "support_generators.R"))
source(file.path(this_dir_, "simulation_utils.R"))


# ==============================================================================
# Monte Carlo base parameters
# ==============================================================================

MC_BASE <- list(
  n         = 200L,
  rho_scale = 1.0,
  snr       = 2.0,
  tFDR      = 0.1,
  K         = 20L,
  num_MC    = 200L,
  seed      = 2026L,
  corr_max  = 0.98,
  hc_dist   = "single"
)


# ==============================================================================
# Part 2: Monte Carlo simulation -- SNR sweep
# ==============================================================================

if (run_part_2) {
  trx_gvs_12_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    .run_snr_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("Quasi-hapgen GVS MC - SNR sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  rho_scale: %.2f\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$rho_scale))
      cat(sprintf("n=%d  p=500  SNR swept  %d MC\n",
                  MC_BASE$n, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(snr_grid, function(snr_val) {
        cat(sprintf("\n  [SNR = %.2f]  running %d MC trials ...\n",
                    snr_val, MC_BASE$num_MC))
        r <- .run_mc_hapgen_snr(
          n         = MC_BASE$n,
          snr       = snr_val,
          rho_scale = MC_BASE$rho_scale,
          tFDR      = MC_BASE$tFDR,
          K         = MC_BASE$K,
          num_MC    = MC_BASE$num_MC,
          seed      = MC_BASE$seed,
          GVS_type  = GVS_type,
          corr_max  = MC_BASE$corr_max,
          hc_dist   = MC_BASE$hc_dist
        )
        c(list(SNR = snr_val), r)
      })
    }

    res_en <- .run_snr_sweep("EN")
    cat("\n  Results - EN:\n")
    .print_table(res_en, "SNR")

    res_ien <- .run_snr_sweep("IEN")
    cat("\n  Results - IEN:\n")
    .print_table(res_ien, "SNR")
  }

  trx_gvs_12_02()

}  # end Part 2
# ==============================================================================


# ==============================================================================
# Part 3: Monte Carlo simulation -- rho_scale sweep
# ==============================================================================

if (run_part_3) {
  trx_gvs_12_03 <- function() {

    rho_scale_grid <- c(0.15, 0.30, 0.50, 0.70, 0.85, 1.00)

    .run_rho_scale_sweep <- function(GVS_type) {
      cat(strrep("=", 70), "\n")
      cat("Quasi-hapgen GVS MC - rho_scale sweep\n")
      cat(sprintf("tFDR: %.2f  |  GVS_type: %s  |  SNR: %.2f\n",
                  MC_BASE$tFDR, GVS_type, MC_BASE$snr))
      cat(sprintf("n=%d  p=500  rho_scale swept  %d MC\n",
                  MC_BASE$n, MC_BASE$num_MC))
      cat(strrep("=", 70), "\n\n")

      lapply(rho_scale_grid, function(rs_val) {
        cat(sprintf("\n  [rho_scale = %.2f]  running %d MC trials ...\n",
                    rs_val, MC_BASE$num_MC))
        r <- .run_mc_hapgen_snr(
          n         = MC_BASE$n,
          snr       = MC_BASE$snr,
          rho_scale = rs_val,
          tFDR      = MC_BASE$tFDR,
          K         = MC_BASE$K,
          num_MC    = MC_BASE$num_MC,
          seed      = MC_BASE$seed,
          GVS_type  = GVS_type,
          corr_max  = MC_BASE$corr_max,
          hc_dist   = MC_BASE$hc_dist
        )
        c(list(rho_sc = rs_val), r)
      })
    }

    res_en <- .run_rho_scale_sweep("EN")
    cat("\n  Results - EN:\n")
    .print_table(res_en, "rho_sc", "%-7.2f")

    res_ien <- .run_rho_scale_sweep("IEN")
    cat("\n  Results - IEN:\n")
    .print_table(res_ien, "rho_sc", "%-7.2f")
  }

  trx_gvs_12_03()

}  # end Part 3
# ==============================================================================


# ==============================================================================
# Part 4: 2D phase-transition sweep: SNR x rho_scale
# ==============================================================================

if (run_part_4) {
  trx_gvs_12_04 <- function() {

    snr_grid       <- c(0.2, 0.5, 1.0, 2.0, 5.0)
    rho_scale_grid <- c(0.30, 0.50, 0.70, 0.85, 1.00)
    n_snr          <- length(snr_grid)
    n_rs           <- length(rho_scale_grid)

    row_nms <- paste0("snr=", sprintf("%.2f", snr_grid))
    col_nms <- paste0("rs=",  sprintf("%.2f", rho_scale_grid))

    mat_TPP_en  <- matrix(NA_real_, nrow = n_snr, ncol = n_rs,
                          dimnames = list(row_nms, col_nms))
    mat_FDP_en  <- mat_TPP_en
    mat_TPP_ien <- mat_TPP_en
    mat_FDP_ien <- mat_TPP_en

    cat(strrep("=", 70), "\n")
    cat("Quasi-hapgen GVS MC  |  Part 4: SNR x rho_scale sweep\n")
    cat(sprintf("n=%d, p=500, %d MC, tFDR=%.2f, corr_max=%.2f\n",
                MC_BASE$n, MC_BASE$num_MC, MC_BASE$tFDR, MC_BASE$corr_max))
    cat(sprintf("  snr_grid       : {%s}\n", paste(snr_grid, collapse = ", ")))
    cat(sprintf("  rho_scale_grid : {%s}\n",
                paste(sprintf("%.2f", rho_scale_grid), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    total_cells <- n_snr * n_rs
    cell_idx    <- 0L

    for (i_snr in seq_len(n_snr)) {
      snr_val <- snr_grid[i_snr]
      for (i_rs in seq_len(n_rs)) {
        rs_val   <- rho_scale_grid[i_rs]
        cell_idx <- cell_idx + 1L
        prefix   <- sprintf("  [%d/%d] SNR=%.2f  rho_scale=%.2f",
                            cell_idx, total_cells, snr_val, rs_val)

        cat(sprintf("%s  [EN]  running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_en <- .run_mc_hapgen_snr(
          n         = MC_BASE$n,
          snr       = snr_val,
          rho_scale = rs_val,
          tFDR      = MC_BASE$tFDR,
          K         = MC_BASE$K,
          num_MC    = MC_BASE$num_MC,
          seed      = MC_BASE$seed,
          GVS_type  = "EN",
          corr_max  = MC_BASE$corr_max,
          hc_dist   = MC_BASE$hc_dist
        )
        mat_TPP_en[i_snr, i_rs] <- r_en$mean_TPP
        mat_FDP_en[i_snr, i_rs] <- r_en$mean_FDP
        cat(sprintf("    %s [EN]  done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_en$mean_TPP, r_en$mean_FDP))

        cat(sprintf("%s  [IEN] running %d MC trials ...\n", prefix, MC_BASE$num_MC))
        r_ien <- .run_mc_hapgen_snr(
          n         = MC_BASE$n,
          snr       = snr_val,
          rho_scale = rs_val,
          tFDR      = MC_BASE$tFDR,
          K         = MC_BASE$K,
          num_MC    = MC_BASE$num_MC,
          seed      = MC_BASE$seed,
          GVS_type  = "IEN",
          corr_max  = MC_BASE$corr_max,
          hc_dist   = MC_BASE$hc_dist
        )
        mat_TPP_ien[i_snr, i_rs] <- r_ien$mean_TPP
        mat_FDP_ien[i_snr, i_rs] <- r_ien$mean_FDP
        cat(sprintf("    %s [IEN] done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r_ien$mean_TPP, r_ien$mean_FDP))
      }
    }

    .print_matrix(mat_TPP_en,  "mean_TPP [EN]  (rows: SNR, cols: rho_scale)")
    .print_matrix(mat_FDP_en,  "mean_FDP [EN]  (rows: SNR, cols: rho_scale)")
    .print_matrix(mat_TPP_ien, "mean_TPP [IEN] (rows: SNR, cols: rho_scale)")
    .print_matrix(mat_FDP_ien, "mean_FDP [IEN] (rows: SNR, cols: rho_scale)")
  }

  trx_gvs_12_04()

}  # end Part 4
# ==============================================================================


cat("\nQuasi-hapgen LD-block GVS MC simulations complete.\n")
