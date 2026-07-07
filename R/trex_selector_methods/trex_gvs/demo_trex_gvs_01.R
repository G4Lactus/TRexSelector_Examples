# =================================================================================
# demo_trex_gvs_01.R
# =================================================================================
#
# T-Rex+GVS single-run demo for the Hastie (Zou & Hastie 2005) DGP.
#
#  Part 1: Single-run demo — data visualization and selector application.
#
# Monte Carlo simulations (Parts 2–4) are in simulation_trex_gvs_01.R.
#
# =================================================================================

library(TRexSelectorNeo)
library(plotly)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))
source(file.path(this_dir_, "dgp_hastie.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

# Part 1 parameters
PARAMS <- list(
  n        = 200L,
  p        = 500L,
  sd_x     = sqrt(0.01),
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

#' Print a formatted single-run result block.
.print_result <- function(scenario_name, dat, result, tFDR) {
  sel_set <- which(result$selected_var == 1L)
  n_sel   <- length(sel_set)
  n_tp    <- length(intersect(sel_set, which(dat$beta != 0)))
  n_fp    <- n_sel - n_tp
  tpp_val <- TRexSelectorNeo::compute_tpp(result$selected_indices, which(dat$beta != 0))
  fdp_val <- TRexSelectorNeo::compute_fdp(result$selected_indices, which(dat$beta != 0))

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
              tpp_val, fdp_val, tFDR))
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  trx_gvs_01_hastie_01 <- function() {

    # Header
    # --------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Hastie GVS Demo  |  Part 1: Single-run\n")
    cat(sprintf("n=%d,  p=%d,  SNR=%.2f,  sd_x=%.4f,  tFDR=%.2f\n",
                PARAMS$n, PARAMS$p, PARAMS$snr, PARAMS$sd_x, PARAMS$tFDR))
    cat(strrep("=", 70), "\n\n")
    # --------------------------------------------

    # Generate data with SNR-calibrated noise
    cat("[Part 1] Generating Hastie data ...\n")
    dat_p1 <- dgp_hastie_snr(
      n    = PARAMS$n,
      p    = PARAMS$p,
      snr  = PARAMS$snr,
      sd_x = PARAMS$sd_x,
      seed = PARAMS$seed
    )
    cat(sprintf("[Part 1] Active variables: %d  |  sigma_y = %.4f\n\n",
                dat_p1$s, dat_p1$sigma_y))

    # Visualize correlation matrix — confirms 3 dense diagonal blocks
    cat("[Part 1] Plotting correlation matrix ...\n")
    plot_cormat(cor(dat_p1$X))

    # Run T-Rex+GVS with EN
    cat("[Part 1] Running trex+GVS (EN) ...\n\n")
    res_en <- TRexSelectorNeo::TRexGVSSelector$new(
      dat_p1$X, dat_p1$y, tFDR = PARAMS$tFDR, seed = -1L, verbose = FALSE,
      gvs_control = TRexSelectorNeo::trex_gvs_control(
        gvs_type = "EN", corr_max = PARAMS$corr_max,
        hc_linkage = cap_hc(PARAMS$hc_dist)),
      control = TRexSelectorNeo::trex_control(solver = "TLARS", K = PARAMS$K))
    res_en$select()
    .print_result("Part 1 \u2014 Hastie GVS  [EN]", dat_p1, res_en, PARAMS$tFDR)

    # Run T-Rex+GVS with IEN
    cat("[Part 1] Running trex+GVS (IEN) ...\n\n")
    res_ien <- TRexSelectorNeo::TRexGVSSelector$new(
      dat_p1$X, dat_p1$y, tFDR = PARAMS$tFDR, seed = -1L, verbose = FALSE,
      gvs_control = TRexSelectorNeo::trex_gvs_control(
        gvs_type = "IEN", corr_max = PARAMS$corr_max,
        hc_linkage = cap_hc(PARAMS$hc_dist)),
      control = TRexSelectorNeo::trex_control(solver = "TLARS", K = PARAMS$K))
    res_ien$select()
    .print_result("Part 1 \u2014 Hastie GVS  [IEN]", dat_p1, res_ien, PARAMS$tFDR)
  }

  trx_gvs_01_hastie_01()

}  # end Part 1
# ==============================================================================


cat("\nHastie GVS demo complete.\n")
