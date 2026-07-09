# ==============================================================================
# demo_trex_da_02_equi.R
# ==============================================================================
# DA-TRex demo and MC simulation for the full-matrix equi-correlated DGP.
#
# Reuses the BT factor-model DGP by setting rho_within == rho_between.
# ==============================================================================

library(TRexSelectorNeo)
library(plotly)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                      error = function(e) {
                        args <- commandArgs(trailingOnly = FALSE)
                        file_arg <- grep("--file=", args, value = TRUE)
                        if (length(file_arg) > 0) {
                          dirname(normalizePath(sub("--file=", "", file_arg[1])))
                        } else {
                          "."
                        }
                      })

# for make_support_random
source(file.path(this_dir_, "..", "..", "support_generators.R"))

# for dgp_bt_snr
source(file.path(this_dir_, "..", "trex_da_dgps.R"))
source(file.path(this_dir_, "..", "..", "simulation_utils.R"))

# ==============================================================================
# Part 1: Single-run demo
# n=150, p=500, s=5, rho=0.25, SNR=2.0
# ==============================================================================
if (FALSE) {

  trx_da_02_equi_01 <- function() {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("  EQUI Demo | Part 1: Single-run\n")
    cat("  n=150 p=500 s=5 rho=0.25 SNR=2.0\n")
    cat(strrep("=", 70), "\n\n")

    support_p1 <- make_support_random(s = 5L, p = 500L, seed = 2026L)

    # Set rho_within and rho_between both to 0.25 to make it pure equi-correlated
    dat_p1 <- dgp_bt_snr(
      n           = 150L,
      p           = 500L,
      support     = support_p1,
      amplitude   = 3.0,
      n_blocks    = 5L,
      rho_within  = 0.25,
      rho_between = 0.25,
      snr         = 2.0,
      seed        = 2026L
    )

    plot_cormat(cor(dat_p1$X))

    res_p1 <- TRexSelectorNeo::TRexDASelector$new(
      dat_p1$X, dat_p1$y, tFDR = 0.2, seed = -1L, verbose = FALSE,
      da_control = TRexSelectorNeo::trex_da_control(da_method = "EQUI"),
      control = TRexSelectorNeo::trex_control(solver = "TLARS", K = 20L))
    res_p1$select()

    cat(sprintf("  Selected variables: %d\n\n", sum(res_p1$selected_var)))

  }
  # Run demo for Part 1
  trx_da_02_equi_01()

} # end of part 1
# ==============================================================================


# ==============================================================================
# Part 2: Monte Carlo simulation — SNR sweep
# Fixed rho_within = 0.25, rho_between = 0.25.
# Swept: SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
# ==============================================================================
if (TRUE) {

  trx_da_02_equi_02 <- function() {
    MC_PARAMS <- list(
      n           = 300L,
      p           = 1000L,
      s           = 10L,
      amplitude   = 3.0,
      rho_within  = 0.25,
      rho_between = 0.25,
      tFDR        = 0.2,
      K           = 20L,
      num_MC      = 200L,
      seed        = 2026L
    )

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    cat(strrep("=", 70), "\n")
    cat(sprintf(" EQUI MC — SNR sweep | rho=%.2f, %d MC\n",
                MC_PARAMS$rho_within, MC_PARAMS$num_MC))
    cat(strrep("=", 70), "\n\n")

    snr_results <- lapply(snr_grid, function(snr_val) {
      cat(sprintf("\n [SNR = %.2f] running %d MC trials ...\n", snr_val, MC_PARAMS$num_MC))
      r <- .run_mc_equi(
        n        = MC_PARAMS$n,
        p        = MC_PARAMS$p,
        support  = NULL,
        amplitude = MC_PARAMS$amplitude,
        s        = MC_PARAMS$s,
        rho      = MC_PARAMS$rho_within,
        n_blocks = 10L,
        snr      = snr_val,
        tFDR     = MC_PARAMS$tFDR,
        K        = MC_PARAMS$K,
        num_MC   = MC_PARAMS$num_MC,
        seed     = MC_PARAMS$seed,
        n_cores  = num_cores
      )
      c(list(SNR = snr_val), r)
    })

    cat("\n\n Results — SNR sweep:\n")
    .print_table(snr_results, "SNR", "%-8.2f")
  }

  # Run demo for Part 2
  trx_da_02_equi_02()

} # end of part 2
# ==============================================================================


# ==============================================================================
# Part 3: Monte Carlo simulation — 1D Pure Equi Sweep
# Fixed SNR = 2.0.
# Swept: rho in seq(0.0, 0.9, by = 0.1).
# ==============================================================================
if (TRUE) {

  trx_da_02_equi_03 <- function() {
    MC_PARAMS_P3 <- list(
      n         = 300L,
      p         = 1000L,
      s         = 10L,
      amplitude = 3.0,
      snr       = 2.0,
      tFDR      = 0.2,
      K         = 20L,
      num_MC    = 200L,
      seed      = 2026L,
      n_blocks  = 10L
    )

    rho_grid <- seq(0.0, 0.9, by = 0.1)

    cat(strrep("=", 70), "\n")
    cat(sprintf(" EQUI MC — 1D Pure Equi Sweep | SNR=%.1f, %d MC\n",
                MC_PARAMS_P3$snr, MC_PARAMS_P3$num_MC))
    cat(" Sweeping rho (where rho_within == rho_between)\n")
    cat(strrep("=", 70), "\n\n")

    rho_results <- lapply(rho_grid, function(rho_val) {
      cat(sprintf("\n [rho = %.1f] running %d MC trials ...\n", rho_val, MC_PARAMS_P3$num_MC))
      r <- .run_mc_equi(
        n        = MC_PARAMS_P3$n,
        p        = MC_PARAMS_P3$p,
        support  = NULL,
        amplitude = MC_PARAMS_P3$amplitude,
        s        = MC_PARAMS_P3$s,
        rho      = rho_val,
        n_blocks = MC_PARAMS_P3$n_blocks,
        snr      = MC_PARAMS_P3$snr,
        tFDR     = MC_PARAMS_P3$tFDR,
        K        = MC_PARAMS_P3$K,
        num_MC   = MC_PARAMS_P3$num_MC,
        seed     = MC_PARAMS_P3$seed,
        n_cores  = num_cores
      )
      c(list(rho = rho_val), r)
    })

    cat("\n\n Results — 1D Pure Equi Sweep:\n")
    .print_table(rho_results, "rho", "%-8.2f")
  }

  # Run demo for Part 3
  trx_da_02_equi_03()

} # end of part 3
# ==============================================================================
