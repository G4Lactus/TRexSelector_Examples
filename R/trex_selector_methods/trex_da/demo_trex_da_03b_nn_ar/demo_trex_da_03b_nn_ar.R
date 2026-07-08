# ==============================================================================
# demo_trex_da_04_nn_ar.R
# ==============================================================================
#
# TREX+DA+NN demo under AR(1) data design.
#
# The data generating process is defined in "dgp_ar1_snr.R".
#
# All p predictors share a banded covariance structure.
#
# Part 1: Single-run demo.
#
# Part 2: Monte Carlo simulation sweeping SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
# Part 3: Monte Carlo simulation sweeping rho in {0.1, 0.2, ..., 0.9}.
#
# Part 4: Monte Carlo simulation sweeping both SNR and rho:
#        rho in {0.1, 0.2, ..., 0.9}
#        SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
#
# ==============================================================================

# Load packages
library(TRexSelectorNeo)
library(plotly)
library(parallel)

num_cores <- 8

# Helper: plot a correlation matrix as a heatmap
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

# Source files
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_ar1_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))

# ------------------------------------------------------------------------------
# NN availability guard
# ------------------------------------------------------------------------------
# trex+DA+NN is implemented in the C++ core (DAMethod::NN) but not yet exposed by
# the installed TRexSelectorNeo R binding (trex_da_control() accepts AR1, EQUI,
# BT only; NN throws "Unknown DAMethod: NN" at select()). Detect this once and
# exit cleanly instead of erroring mid-run. Remove once the binding exposes NN.
.nn_supported <- tryCatch({
  .x_ <- matrix(stats::rnorm(20L * 8L), 20L, 8L)
  .s_ <- TRexSelectorNeo::TRexDASelector$new(
    .x_, as.numeric(.x_[, 1L]), tFDR = 0.5, seed = 1L, verbose = FALSE,
    da_control = TRexSelectorNeo::trex_da_control(da_method = "NN"),
    control    = TRexSelectorNeo::trex_control(K = 5L))
  .s_$select()
  TRUE
}, error = function(e) FALSE)

if (!.nn_supported) {
  cat("\n[SKIP] trex+DA+NN is not exposed by the installed TRexSelectorNeo ",
      "R binding.\n",
      "       The C++ core supports DAMethod::NN; this demo will run once the\n",
      "       R binding adds it. Skipping (nothing to compute).\n\n", sep = "")
  if (!interactive()) quit(save = "no", status = 0)
}

# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n         = 150L,   # observations (Part 1)
  p         = 500L,   # predictors (Part 1)
  s         = 5L,     # active predictors
  amplitude = 3.0,    # signal coefficient for each active beta_j
  snr       = 2.0,    # signal-to-noise ratio
  kappa     = 5L,     # MA bandwidth
  rho       = 0.7,    # MA decay parameter
  K         = 20L,    # random experiments per selector run
  tFDR      = 0.2,    # target FDR level
  seed      = 2026L,  # global RNG seed
  max_gap   = 20L     # max gap for CappedSpread support
)

MC_BASE <- list(
  n         = 300L,   # observations (Parts 2-4)
  p         = 1000L,  # predictors (Parts 2-4)
  s         = 10L,    # active predictors
  amplitude = 3.0,
  kappa     = 3L,
  rho       = 0.7,
  snr       = 2.0,
  tFDR      = 0.2,
  K         = 20L,
  num_MC    = 200L,
  seed      = 2026L,
  max_gap   = 20L
)


# ==============================================================================
# Helpers
# ==============================================================================

#' Print a formatted single-run result block.
.print_result <- function(scenario_name, dat, result, tFDR) {
  sel_set <- which(result$selected_var == 1L)
  n_sel   <- length(sel_set)
  n_tp    <- length(intersect(sel_set, dat$true_support))
  n_fp    <- n_sel - n_tp
  tpp_val <- TRexSelectorNeo::compute_tpp(result$selected_indices, dat$true_support)
  fdp_val <- TRexSelectorNeo::compute_fdp(result$selected_indices, dat$true_support)

  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Data: n=%d  p=%d  s=%d  tFDR=%.2f\n",
              dat$n, dat$p, dat$s, tFDR))
  cat(sprintf("        amplitude=%.2f  SNR=%.2f  kappa=%d  rho=%.2f\n",
              PARAMS$amplitude, dat$snr, dat$kappa, dat$rho))
  cat(sprintf("  True support (1-based): {%s}\n",
              paste(dat$true_support, collapse = ", ")))
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
# n=150, p=500, s=5, rho=0.7, SNR=2.0. CappedSpread (max_gap=20).
# ==============================================================================

if (FALSE) {

  trx_da_04_ar_01 <- function() {

    # Generate support
    # support_p1 <- make_support_capped_spread(PARAMS$s, PARAMS$p, PARAMS$max_gap)
    support_p1 <- make_support_random(PARAMS$s, PARAMS$p, sample(150, 1))

    # Header
    # -------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf(" AR(1) Demo | Part 1: Single-run\n"))
    cat(sprintf(
      " n=%d p=%d s=%d rho=%.1f SNR=%.1f\n",
      PARAMS$n, PARAMS$p, PARAMS$s, PARAMS$rho, PARAMS$snr
    ))
    cat(sprintf(" support: Random -> {%s}\n",
                paste(support_p1, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")
    # -------------------------------------------------

    # Data generation
    cat("[Part 1] Generating AR(1) data ...\n")
    dat_p1 <- dgp_ar1_snr(
      n = PARAMS$n,
      p = PARAMS$p,
      support = support_p1,
      amplitude = PARAMS$amplitude,
      rho = PARAMS$rho,
      snr = PARAMS$snr,
      seed = NULL # PARAMS$seed
    )

    # Correlation structure of X
    plot_cormat(cor(dat_p1$X))

    # Run TRex+DA+NN
    cat("[Part 1] Running trex+DA+NN ...\n\n")
    res_p1 <- TRexSelectorNeo::TRexDASelector$new(
      dat_p1$X, dat_p1$y, tFDR = PARAMS$tFDR, seed = -1L, verbose = FALSE,
      da_control = TRexSelectorNeo::trex_da_control(da_method = "NN"),
      control = TRexSelectorNeo::trex_control(solver = "TLARS", K = PARAMS$K))
    res_p1$select()

    # Print results
    .print_result(
      "Part 1 — AR(1) Demo [trex+DA+NN, rho=0.7]",
      dat_p1, res_p1, PARAMS$tFDR
    )
  }

  # Run demo for Part 1
  trx_da_04_ar_01()

} # end Part 1
# ==============================================================================




# ==============================================================================
# Part 2: MC simulation — SNR sweep
# Fixed: n=300, p=1000, s=10, rho=0.7, tFDR=0.2, K=20, 200 MC.
# Swept: SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
# Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================

if (FALSE) {

  trx_da_04_ar_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    .run_snr_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      # Header
      # -------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("AR(1) MC — SNR sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d p=%d s=%d rho=%.1f SNR swept %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$rho, MC_BASE$num_MC))
      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # -------------------------------------------------

      lapply(snr_grid, function(snr_val) {
        cat(sprintf("\n [SNR = %.2f] running %d MC trials ...\n",
                    snr_val, MC_BASE$num_MC))
        r <- .run_mc(
          n              = MC_BASE$n,
          p              = MC_BASE$p,
          support        = support,
          amplitude      = MC_BASE$amplitude,
          snr            = snr_val,
          tFDR           = MC_BASE$tFDR,
          K              = MC_BASE$K,
          num_MC         = MC_BASE$num_MC,
          seed           = MC_BASE$seed,
          dgp_fun        = dgp_ar1_snr,
          dgp_extra_args = list(rho = MC_BASE$rho),
          trex_method    = "trex+DA+NN"
        )
        c(list(SNR = snr_val), r)
      })
    }

    res_capped <- .run_snr_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_table(res_capped, "SNR")

    res_random <- .run_snr_sweep("random", "Random (per-trial)")
    cat("\n Results — Random support (redrawn per trial):\n")
    .print_table(res_random, "SNR")
  }

  trx_da_04_ar_02()

} # end Part 2
# ==============================================================================


# ==============================================================================
# Part 3: MC simulation — Rho sweep
# Fixed: n=300, p=1000, s=10, SNR=2.0, tFDR=0.2, K=20, 200 MC.
# Swept: rho in {0.1, 0.2, ..., 0.9}.
#
# Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================

if (FALSE) {

  trx_da_04_ar_03 <- function() {

    rho_grid <- seq(0.1, 0.9, by = 0.1)

    .run_rho_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      cat(strrep("=", 70), "\n")
      cat("AR(1) MC — Rho sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d p=%d s=%d rho swept SNR=%.1f %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$snr, MC_BASE$num_MC))
      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")

      lapply(rho_grid, function(rho_val) {
        cat(sprintf("\n [rho = %.1f] running %d MC trials ...\n",
                    rho_val, MC_BASE$num_MC))
        r <- .run_mc(
          n              = MC_BASE$n,
          p              = MC_BASE$p,
          support        = support,
          amplitude      = MC_BASE$amplitude,
          snr            = MC_BASE$snr,
          tFDR           = MC_BASE$tFDR,
          K              = MC_BASE$K,
          num_MC         = MC_BASE$num_MC,
          seed           = MC_BASE$seed,
          dgp_fun        = dgp_ar1_snr,
          dgp_extra_args = list(rho = rho_val),
          trex_method    = "trex+DA+NN"
        )
        c(list(rho = rho_val), r)
      })
    }

    res_capped <- .run_rho_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_table(res_capped, "rho", "%-6.1f")

    res_random <- .run_rho_sweep("random", "Random (per-trial)")
    cat("\n Results — Random support (redrawn per trial):\n")
    .print_table(res_random, "rho", "%-6.1f")
  }

  trx_da_04_ar_03()

} # end Part 3
# ==============================================================================


# ==============================================================================
# Part 4: Monte Carlo simulation — 2D SNR x Rho sweep
# Fixed: n=300, p=1000, s=10, tFDR=0.2, K=20, 200 MC.
# Swept: rho in {0.1, 0.2, ..., 0.9}
#        SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
#
# Goal: identify the narrow operating region of DA+NN under AR(1).
# Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================

if (TRUE) {

  trx_da_04_ar_04 <- function() {

    snr_grid <- c(2.0, 5.0) #c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)
    rho_grid <- seq(0.1, 0.9, by = 0.1)

    .run_snr_rho_sweep <- function(support_type, support_label) {

      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      n_s <- length(snr_grid)
      n_r <- length(rho_grid)

      row_nms <- paste0("SNR=", sprintf("%.2f", snr_grid))
      col_nms <- paste0("rho=", sprintf("%.1f", rho_grid))

      mat_TPP <- matrix(NA_real_, nrow = n_s, ncol = n_r,
                        dimnames = list(row_nms, col_nms))
      mat_FDP <- matrix(NA_real_, nrow = n_s, ncol = n_r,
                        dimnames = list(row_nms, col_nms))

      cat(strrep("=", 70), "\n")
      cat("AR(1) MC — 2D SNR x Rho sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | %d MC\n",
                  MC_BASE$tFDR, support_label, MC_BASE$num_MC))
      cat(" Sweeping SNR (rows) and rho (cols)\n")
      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")

      total_cells <- n_s * n_r
      cell_idx <- 0L

      for (i_s in seq_len(n_s)) {
        snr_val <- snr_grid[i_s]

        for (i_r in seq_len(n_r)) {
          rho_val <- rho_grid[i_r]
          cell_idx <- cell_idx + 1L

          prefix <- sprintf(" [%2d/%2d] SNR=%-4.2f rho=%.1f",
                            cell_idx, total_cells, snr_val, rho_val)
          cat(sprintf("%s running %d MC trials ...\n", prefix, MC_BASE$num_MC))

          r <- .run_mc(
            n              = MC_BASE$n,
            p              = MC_BASE$p,
            support        = support,
            amplitude      = MC_BASE$amplitude,
            snr            = snr_val,
            tFDR           = MC_BASE$tFDR,
            K              = MC_BASE$K,
            num_MC         = MC_BASE$num_MC,
            seed           = MC_BASE$seed,
            dgp_fun        = dgp_ar1_snr,
            dgp_extra_args = list(rho = rho_val),
            trex_method    = "trex+DA+NN"
          )

          mat_TPP[i_s, i_r] <- r$mean_TPP
          mat_FDP[i_s, i_r] <- r$mean_FDP

          cat(sprintf(" %s done. TPP=%.3f FDP=%.3f\n\n",
                      prefix, mat_TPP[i_s, i_r], mat_FDP[i_s, i_r]))
        }
      }

      list(mean_TPP = mat_TPP, mean_FDP = mat_FDP)
    }

    # --- Run A: CappedSpread support (fixed, max_gap=20) ---
    if (FALSE) {
      res_capped <- .run_snr_rho_sweep(
        "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
      )
      cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
      .print_matrix(res_capped$mean_TPP, "mean_TPP (rows: SNR, cols: rho)")
      .print_matrix(res_capped$mean_FDP, "mean_FDP (rows: SNR, cols: rho)")
    }

    # --- Run B: Random support (redrawn per trial) ---
    res_random <- .run_snr_rho_sweep("random", "Random (per-trial)")
    cat("\n Results — Random support (redrawn per trial):\n")
    .print_matrix(res_random$mean_TPP, "mean_TPP (rows: SNR, cols: rho)")
    .print_matrix(res_random$mean_FDP, "mean_FDP (rows: SNR, cols: rho)")

    invisible(list(
      capped = res_capped,
      random = res_random
    ))
  }

  # Run demo for Part 4
  trx_da_04_ar_04()

} # end Part 4
# ==============================================================================
