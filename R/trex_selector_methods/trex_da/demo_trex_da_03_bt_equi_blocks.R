# ==============================================================================
# demo_trex_da_03_bt_equi_blocks.R
# ==============================================================================
#
# DA-TRex demo and MC simulation for the Binary-Tree (BT) hierarchical-block data
# model with equi-correlated blocks.
#
#  Part 1: Single-run demo.
#
#  Part 2: MC simulation — SNR sweep, factor-model DGP, average linkage.
#
#  Part 3: MC simulation — SNR sweep, factor-model DGP, linkage sweep over
#          hc_dist in {"single", "complete", "average"}.
#
# ==============================================================================

# Resolve the directory of this file, works for both source() and Rscript.
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_bt_snr.R"))

library(TRexSelectorNeo)
library(plotly)
library(parallel)

num_cores <- 6L

source(file.path(this_dir_, "..", "simulation_utils.R"))


# ==============================================================================
#  Part 1: Single-run demo
# ==============================================================================
#
# Evaluation helpers
# ==============================================================================

# Print a formatted single-run result block.
.print_result <- function(scenario_name, dat, result, tFDR) {
  sel_set <- which(result$selected_var == 1L)
  n_sel   <- length(sel_set)
  n_tp    <- length(intersect(sel_set, dat$true_support))
  n_fp    <- n_sel - n_tp
  tpp_val <- TRexSelectorNeo::compute_tpp(result$selected_indices, dat$true_support)
  fdp_val <- TRexSelectorNeo::compute_fdp(result$selected_indices, dat$true_support)

  # Header
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf(
    "  Data:       n = %d,  p = %d,  s = %d,  SNR = %.2f\n",
    dat$n, dat$p, dat$s, dat$snr
  ))
  cat(sprintf(
    "  BT:         %d blocks,  rho_within = %.2f,  rho_between = %.2f\n",
    dat$n_blocks, dat$rho_within, dat$rho_between
  ))
  cat(sprintf(
    "  True support (1-based): {%s}\n",
    paste(dat$true_support, collapse = ", ")
  ))
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

#  Part 1: Single-run demo
# ==============================================================================
if (FALSE) {
  trx_da_03_bt_equi_blks_01 <- function() {
    PARAMS <- list(
      n           = 150L,   # observations
      p           = 500L,   # predictors
      s           = 5L,     # active predictors (= n_blocks: one active per block)
      amplitude   = 3.0,    # signal coefficient for each active beta_j
      snr         = 2.0,    # signal-to-noise ratio
      n_blocks    = 5L,     # number of hierarchical blocks (must equal s)
      rho_within  = 0.2,    # within-block correlation
      rho_between = 0.05,   # between-block correlation
      K           = 20L,    # random experiments per selector run
      tFDR        = 0.2,    # target FDR level
      seed        = 2026L   # global RNG seed (= base_seed in C++)
    )

    # Header
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf(
      "  BT Demo  |  n=%d  p=%d  s=%d  blocks=%d  rho_w=%.1f  rho_b=%.1f  support=OnePerBlock\n",
      PARAMS$n, PARAMS$p, PARAMS$s, PARAMS$n_blocks,
      PARAMS$rho_within, PARAMS$rho_between
    ))
    cat(strrep("=", 70), "\n\n")

    # One-per-block support: first predictor of each block.
    # block_size = 500/5 = 100 -> support = {1, 101, 201, 301, 401}
    support_p1 <- make_support_bt_one_per_block(
      p        = PARAMS$p,
      n_blocks = PARAMS$n_blocks
    )

    cat("[Part 1] Generating BT data\n")
    cat(sprintf("  n=%d, p=%d, s=%d, blocks=%d, rho_w=%.1f, rho_b=%.1f, SNR=%.1f ...\n",
                PARAMS$n, PARAMS$p, PARAMS$s, PARAMS$n_blocks,
                PARAMS$rho_within, PARAMS$rho_between, PARAMS$snr))

    dat_p1 <- dgp_bt_snr(
      n           = PARAMS$n,
      p           = PARAMS$p,
      support     = support_p1,
      amplitude   = PARAMS$amplitude,
      n_blocks    = PARAMS$n_blocks,
      rho_within  = PARAMS$rho_within,
      rho_between = PARAMS$rho_between,
      snr         = PARAMS$snr,
      seed        = PARAMS$seed
    )

    cor_mat_d_BT <- cor(dat_p1$X)
    plot_cormat(cor_mat_d_BT)

    # Run T-Rex+DA+BT
    cat("[Part 1] Running trex+DA+BT (hc_dist = 'single') ...\n\n")
    res_p1 <- TRexSelectorNeo::TRexDASelector$new(
      dat_p1$X, dat_p1$y, tFDR = PARAMS$tFDR, seed = -1L, verbose = FALSE,
      da_control = TRexSelectorNeo::trex_da_control(
        da_method = "BT", hc_linkage = cap_hc("single")),
      control = TRexSelectorNeo::trex_control(solver = "TLARS", K = PARAMS$K))
    res_p1$select()

    .print_result(
      "Part 1 — BT Demo  [trex+DA+BT, TLARS, hc_dist=single]",
      dat_p1, res_p1, PARAMS$tFDR
    )
  }

  # Run demo for Part 1
  trx_da_03_bt_equi_blks_01()
}
# ==============================================================================


# ==============================================================================
#  Part 2: Monte Carlo simulation — SNR sweep
#
#  n=300, p=1000, s=10, n_blocks=10 (one active per block), 200 MC.
#  Support: OnePerBlock (first predictor of each block) -> {1, 101, ..., 901}.
#  SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
# ==============================================================================
if (TRUE) {
  trx_da_03_bt_equi_blks_02 <- function() {
    MC_PARAMS <- list(
      n           = 300L,
      p           = 1000L,
      s           = 10L,     # must equal n_blocks
      amplitude   = 3.0,
      n_blocks    = 10L,
      rho_within  = 0.5,
      rho_between = 0.1,
      tFDR        = 0.2,
      K           = 20L,
      num_MC      = 200L,
      seed        = 2026L
    )

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    # One-per-block support: block_size = 1000/10 = 100 -> {1, 101, 201, ..., 901}
    support_mc <- make_support_bt_one_per_block(
      p        = MC_PARAMS$p,
      n_blocks = MC_PARAMS$n_blocks
    )

    # Header
    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  BT MC — SNR sweep  |  n=%d  p=%d  s=%d  blocks=%d  rho_w=%.1f  rho_b=%.1f  %d MC\n",
      MC_PARAMS$n, MC_PARAMS$p, MC_PARAMS$s, MC_PARAMS$n_blocks,
      MC_PARAMS$rho_within, MC_PARAMS$rho_between, MC_PARAMS$num_MC
    ))
    cat(sprintf("  Support: OnePerBlock  -> {%s, ...}\n",
                paste(head(support_mc, 4), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    snr_results <- lapply(snr_grid, function(snr_val) {
      cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, MC_PARAMS$num_MC))
      r <- .run_mc_bt(
        n           = MC_PARAMS$n,
        p           = MC_PARAMS$p,
        support     = support_mc,
        amplitude   = MC_PARAMS$amplitude,
        n_blocks    = MC_PARAMS$n_blocks,
        rho_within  = MC_PARAMS$rho_within,
        rho_between = MC_PARAMS$rho_between,
        snr         = snr_val,
        tFDR        = MC_PARAMS$tFDR,
        K           = MC_PARAMS$K,
        num_MC      = MC_PARAMS$num_MC,
        seed        = MC_PARAMS$seed,
        hc_dist     = "single",
        n_cores     = num_cores
      )
      c(list(snr = snr_val), r)
    })

    cat("\n")
    .print_table(snr_results, "snr", "%-8.2f")
  }

  # Run demo for Part 2
  trx_da_03_bt_equi_blks_02()
}
# ==============================================================================


# ==============================================================================
#  Part 3: Monte Carlo simulation — linkage sweep over factor-model DGP
#
#  Outer loop sweeps hc_dist over {"single", "complete", "average"} and prints
#  one table per linkage method.
# ==============================================================================
if (TRUE) {
  trx_da_03_bt_equi_blks_03 <- function() {
    MC3_PARAMS <- list(
      n           = 300L,
      p           = 1000L,
      s           = 10L,
      amplitude   = 3.0,
      n_blocks    = 10L,
      rho_within  = 0.8,
      rho_between = 0.1,
      tFDR        = 0.2,
      K           = 20L,
      num_MC      = 200L,
      seed        = 2026L
    )

    snr_grid3    <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
    linkage_grid <- c("single", "complete", "average")

    support_mc3 <- make_support_bt_one_per_block(
      p        = MC3_PARAMS$p,
      n_blocks = MC3_PARAMS$n_blocks
    )

    # Header
    # -------------------------------------------------
    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  BT MC — Linkage sweep  |  n=%d  p=%d  s=%d  blocks=%d  rho_w=%.1f  rho_b=%.1f  %d MC\n",
      MC3_PARAMS$n, MC3_PARAMS$p, MC3_PARAMS$s, MC3_PARAMS$n_blocks,
      MC3_PARAMS$rho_within, MC3_PARAMS$rho_between, MC3_PARAMS$num_MC
    ))
    cat(sprintf("  Support: OnePerBlock  -> {%s, ...}\n",
                paste(head(support_mc3, 4), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")
    # -------------------------------------------------

    # Loop over linkage methods, then SNR values, and run MC simulations
    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      lnk_results <- lapply(snr_grid3, function(snr_val) {
        cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, MC3_PARAMS$num_MC))
        r <- .run_mc_bt(
          n           = MC3_PARAMS$n,
          p           = MC3_PARAMS$p,
          support     = support_mc3,
          amplitude   = MC3_PARAMS$amplitude,
          n_blocks    = MC3_PARAMS$n_blocks,
          rho_within  = MC3_PARAMS$rho_within,
          rho_between = MC3_PARAMS$rho_between,
          snr         = snr_val,
          tFDR        = MC3_PARAMS$tFDR,
          K           = MC3_PARAMS$K,
          num_MC      = MC3_PARAMS$num_MC,
          seed        = MC3_PARAMS$seed,
          hc_dist     = lnk,
          n_cores     = num_cores
        )
        c(list(snr = snr_val), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table(lnk_results, "snr", "%-8.2f")
    }
  }

  # Run demo for Part 3
  trx_da_03_bt_equi_blks_03()
}
# ==============================================================================
