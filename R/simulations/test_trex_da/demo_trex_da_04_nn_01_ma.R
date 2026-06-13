# ==============================================================================
# demo_trex_da_04_nn.R
# ==============================================================================
#
# DA-TRex demo and MC simulation for the Nearest-Neighbours (NN) / MA(kappa) DGP.
#
# The data generating process is defined in "dgp_nn_snr.R".
#
# All p predictors share a single MA(kappa) banded covariance structure:
# Sigma[j, k] != 0 only for |j - k| <= kappa.
#
#  Part 1: Single-run demo.
#
#  Part 2: SNR sweep of SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#          CappedSpread vs Random support.
#
#  Part 3: Rho sweep in {0.0, 0.1, ..., 0.9}.
#          CappedSpread vs Random support.
#
#  Part 4: Kappa sweep in {1, 2, 3, 5, 7, 10, 15}.
#          Fixed rho=0.7, SNR=2. CappedSpread support.
#          Explores how the bandwidth of the correlation structure affects
#          the DA+NN selector.
#
#  Part 5: 2D Kappa x Rho sweep. Fixed SNR=2.0.
#          kappa in {1, 2, 3, 5, 7, 10, 15}
#          rho in seq(0.1, 0.9, by = 0.2).
#
# Notes:
# Set rho high between 0.8 and 0.95. This mimics strong sector correlations in financial
# markets. Active variables are most realistically set by CappedSpread, which simulates
# that one true signal is active, e.g., the sector leader.
# If actives are too densely placed in the same kappa-band, they become unidentifiable.
#
# ==============================================================================

# Load packages
library(TRexSelector)
library(plotly)
library(parallel)

num_cores <- 6

# Helper: plot a correlation matrix as a heatmap
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) "."
)

# Source files
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_ar1_snr.R"))
source(file.path(this_dir_, "dgp_nn_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))

# ==============================================================================
# Base parameters
# ==============================================================================

PARAMS <- list(
  n         = 150L,   # observations
  p         = 500L,   # predictors
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
  n         = 300L,   # observations
  p         = 1000L,  # predictors
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
  tpp_val <- TRexSelector::TPP(result$selected_var, dat$beta)
  fdp_val <- TRexSelector::FDP(result$selected_var, dat$beta)

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
  cat(sprintf("  Calibration:  T_stop = %d,  dummies = %d\n",
              result$T_stop, result$num_dummies))
  cat(sprintf("  Selection:    %d selected  |  TP = %d  FP = %d\n",
              n_sel, n_tp, n_fp))
  cat(sprintf("  Rates:        TPP = %.3f  |  FDP = %.3f  (target tFDR <= %.2f)\n",
              tpp_val, fdp_val, tFDR))
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)

}




# ==============================================================================
#  Part 1: Single-run demo
#  n=150, p=500, s=5, kappa=3, rho=0.7, SNR=2.0. CappedSpread (max_gap=20).
# ==============================================================================
if (TRUE) {

  trx_da_04_nn_01 <- function() {

    # Generate support
    # support_p1 <- make_support_capped_spread(PARAMS$s, PARAMS$p, PARAMS$max_gap)
    support_p1 <- make_support_random(PARAMS$s, PARAMS$p, sample(150, 1))


    # Header
    # -------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf("  NN Demo  |  Part 1: Single-run\n"))
    cat(sprintf(
      "  n=%d  p=%d  s=%d  kappa=%d  rho=%.1f  SNR=%.1f\n",
      PARAMS$n, PARAMS$p, PARAMS$s, PARAMS$kappa, PARAMS$rho, PARAMS$snr
    ))
    cat(sprintf("  support: Random  -> {%s}\n",
                paste(support_p1, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")
    # -------------------------------------------------

    # Data generation
    cat("[Part 1] Generating MA(kappa) data ...\n")
    dat_p1 <- dgp_nn_snr(
      n         = PARAMS$n,
      p         = PARAMS$p,
      support   = support_p1,
      amplitude = PARAMS$amplitude,
      kappa     = 1, #PARAMS$kappa,
      rho       = PARAMS$rho,
      snr       = PARAMS$snr,
      seed      = NULL #PARAMS$seed
    )

    # Correlation structure of X
    cor_mat_p1 <- cor(dat_p1$X)
    plot_cormat(cor_mat_p1)


    # Run TRex+DA+NN
    cat("[Part 1] Running trex+DA+NN ...\n\n")
    res_p1 <- TRexSelector::trex(
      X       = dat_p1$X,
      y       = dat_p1$y,
      tFDR    = PARAMS$tFDR,
      K       = PARAMS$K,
      method  = "trex+DA+NN",
      verbose = FALSE
    )

    # Print results
    .print_result(
      "Part 1 — NN Demo  [trex+DA+NN, kappa=3, rho=0.7]",
      dat_p1, res_p1, PARAMS$tFDR
    )
  }

  # Run demo for Part 1
  trx_da_04_nn_01()

}  # end Part 1
# ==============================================================================


# ==============================================================================
#  Part 2: MC simulation — SNR sweep
#  Fixed: n=300, p=1000, s=10, kappa=3, rho=0.7, tFDR=0.2, K=20, 200 MC.
#  Swept: SNR in {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
#
#  Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================
if (FALSE) {

  trx_da_04_nn_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    # Helper: run SNR sweep for one support type.
    .run_snr_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL  # redrawn per trial inside .run_mc_nn()
      }

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("NN MC — SNR sweep\n")
      cat(sprintf("tFDR: %.2f  |  Support: %s  |  amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d  p=%d  s=%d  kappa=%d  rho=%.1f  SNR swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$kappa,
                  MC_BASE$rho, MC_BASE$num_MC))
      if (!is.null(support)) {
        cat(sprintf("  Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat("  Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # ------------------------------------------------

      # Run MC trials for each SNR value.
      lapply(snr_grid, function(snr_val) {
        cat(sprintf("\n  [SNR = %.2f]  running %d MC trials ...\n",
                    snr_val, MC_BASE$num_MC))
        r <- .run_mc_nn(MC_BASE$n, MC_BASE$p, support, MC_BASE$amplitude,
                        MC_BASE$kappa, MC_BASE$rho, snr_val,
                        MC_BASE$tFDR, MC_BASE$K, MC_BASE$num_MC, MC_BASE$seed,
                        n_cores = num_cores)
        c(list(SNR = snr_val), r)
      })
    }

    # --- Run A: CappedSpread support (fixed, max_gap=20) ---
    res_capped <- .run_snr_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n  Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_table(res_capped, "SNR")

    # --- Run B: Random support (redrawn per trial) ---
    res_random <- .run_snr_sweep("random", "Random (per-trial)")
    cat("\n  Results — Random support (redrawn per trial):\n")
    .print_table(res_random, "SNR")
  }

  # Run demo for Part 2
  trx_da_04_nn_02()

}  # end Part 2
# ==============================================================================


# ==============================================================================
#  Part 3: MC simulation — Rho sweep
#  Fixed: n=300, p=1000, s=10, kappa=3, SNR=2.0, tFDR=0.2, K=20, 200 MC.
#  Swept: rho in {0.1, 0.2, ..., 0.9}.
#  Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================
if (FALSE) {

  trx_da_04_nn_03 <- function() {

    rho_grid <- seq(0.1, 0.9, by = 0.1)

    # Helper: run rho sweep for one support type.
    .run_rho_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("NN MC — Rho sweep\n")
      cat(sprintf("tFDR: %.2f  |  Support: %s  |  amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d  p=%d  s=%d  kappa=%d  rho swept  SNR=%.1f  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$kappa,
                  MC_BASE$snr, MC_BASE$num_MC))
      if (!is.null(support)) {
        cat(sprintf("  Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat("  Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # ------------------------------------------------

      # Run MC trials for each rho value.
      lapply(rho_grid, function(rho_val) {
        cat(sprintf("\n  [rho = %.1f]  running %d MC trials ...\n",
                    rho_val, MC_BASE$num_MC))
        r <- .run_mc_nn(MC_BASE$n, MC_BASE$p, support, MC_BASE$amplitude,
                        MC_BASE$kappa, rho_val, MC_BASE$snr,
                        MC_BASE$tFDR, MC_BASE$K, MC_BASE$num_MC, MC_BASE$seed,
                        n_cores = num_cores)
        c(list(rho = rho_val), r)
      })
    }

    # --- Run A: CappedSpread support (fixed, max_gap=20) ---
    res_capped <- .run_rho_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n  Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_table(res_capped, "rho", "%-6.1f")

    # --- Run B: Random support (redrawn per trial) ---
    res_random <- .run_rho_sweep("random", "Random (per-trial)")
    cat("\n  Results — Random support (redrawn per trial):\n")
    .print_table(res_random, "rho", "%-6.1f")
  }

  # Run demo for Part 3
  trx_da_04_nn_03()

}  # end Part 3
# ==============================================================================


# ==============================================================================
#  Part 4: MC simulation — Kappa sweep
#
#  Fixed: n=300, p=1000, s=10, rho=0.7, SNR=2.0, tFDR=0.2, K=20, 200 MC.
#  Swept: kappa in {1, 2, 3, 5, 7, 10, 15}.
#
#  Benchmarking Support policies
#  CappedSpread support (max_gap=20).
#  Random support (redrawn per trial).
#
#  As kappa increases, the banding of the correlation structure widens and
#  the DA+NN method sees more correlated neighbours per variable.
# ==============================================================================
if (FALSE) {

  trx_da_04_nn_04 <- function() {

    kappa_grid <- c(1L, 2L, 3L, 5L, 7L, 10L, 15L)

    .run_kappa_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("NN MC — Kappa sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d p=%d s=%d rho=%.1f SNR=%.1f kappa swept %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$rho,
                  MC_BASE$snr, MC_BASE$num_MC))
      cat(sprintf(" kappa in {%s}\n", paste(kappa_grid, collapse = ", ")))

      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # ------------------------------------------------

      lapply(kappa_grid, function(kappa_val) {
        cat(sprintf("\n [kappa = %d] running %d MC trials ...\n",
                    kappa_val, MC_BASE$num_MC))
        r <- .run_mc_nn(MC_BASE$n, MC_BASE$p, support, MC_BASE$amplitude,
                        kappa_val, MC_BASE$rho, MC_BASE$snr,
                        MC_BASE$tFDR, MC_BASE$K, MC_BASE$num_MC, MC_BASE$seed,
                        n_cores = num_cores)
        c(list(kappa = kappa_val), r)
      })
    }

    res_capped <- .run_kappa_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_table(res_capped, "kappa", "%-6d")

    res_random <- .run_kappa_sweep("random", "Random (per-trial)")
    cat("\n Results — Random support (redrawn per trial):\n")
    .print_table(res_random, "kappa", "%-6d")
  }

  # Run demo for Part 4
  trx_da_04_nn_04()

} # end Part 4
# ==============================================================================


# ==============================================================================
# Part 5: Monte Carlo simulation — 2D Kappa x Rho sweep
# Fixed: n=300, p=1000, s=10, SNR=2.0, tFDR=0.2, K=20, 200 MC.
# Swept: kappa in {1, 2, 3, 5, 7, 10, 15}
#        rho   in {0.1, 0.3, 0.5, 0.7, 0.8, 0.9}
#
# Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================

if (FALSE) {

  trx_da_04_nn_05 <- function() {

    kappa_grid <- c(1L, 2L, 3L, 5L, 7L, 10L, 15L)
    rho_grid   <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9)

    .run_kappa_rho_sweep <- function(support_type, support_label) {

      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      n_k <- length(kappa_grid)
      n_r <- length(rho_grid)

      row_nms <- paste0("kappa=", kappa_grid)
      col_nms <- paste0("rho=", sprintf("%.1f", rho_grid))

      mat_TPP <- matrix(NA_real_, nrow = n_k, ncol = n_r,
                        dimnames = list(row_nms, col_nms))
      mat_FDP <- matrix(NA_real_, nrow = n_k, ncol = n_r,
                        dimnames = list(row_nms, col_nms))

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("NN MC — 2D Kappa x Rho sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | SNR=%.1f | %d MC\n",
                  MC_BASE$tFDR, support_label, MC_BASE$snr, MC_BASE$num_MC))
      cat(" Sweeping kappa (rows) and rho (cols)\n")
      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # ------------------------------------------------

      total_cells <- n_k * n_r
      cell_idx <- 0L

      for (i_k in seq_len(n_k)) {
        kappa_val <- kappa_grid[i_k]

        for (i_r in seq_len(n_r)) {
          rho_val <- rho_grid[i_r]
          cell_idx <- cell_idx + 1L

          prefix <- sprintf(" [%2d/%2d] kappa=%-2d rho=%.1f",
                            cell_idx, total_cells, kappa_val, rho_val)
          cat(sprintf("%s running %d MC trials ...\n", prefix, MC_BASE$num_MC))

          r <- .run_mc_nn(
            n = MC_BASE$n,
            p = MC_BASE$p,
            support = support,
            amplitude = MC_BASE$amplitude,
            kappa = kappa_val,
            rho = rho_val,
            snr = MC_BASE$snr,
            tFDR = MC_BASE$tFDR,
            K = MC_BASE$K,
            num_MC = MC_BASE$num_MC,
            seed = MC_BASE$seed,
            n_cores = num_cores
          )

          mat_TPP[i_k, i_r] <- r$mean_TPP
          mat_FDP[i_k, i_r] <- r$mean_FDP

          cat(sprintf(" %s done. TPP=%.3f FDP=%.3f\n\n",
                      prefix, mat_TPP[i_k, i_r], mat_FDP[i_k, i_r]))
        }
      }

      list(mean_TPP = mat_TPP, mean_FDP = mat_FDP)
    }

    # --- Run A: CappedSpread support (fixed, max_gap=20) ---
    res_capped <- .run_kappa_rho_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_matrix(res_capped$mean_TPP, "mean_TPP (rows: kappa, cols: rho)")
    .print_matrix(res_capped$mean_FDP, "mean_FDP (rows: kappa, cols: rho)")

    # --- Run B: Random support (redrawn per trial) ---
    res_random <- .run_kappa_rho_sweep("random", "Random (per-trial)")
    cat("\n Results — Random support (redrawn per trial):\n")
    .print_matrix(res_random$mean_TPP, "mean_TPP (rows: kappa, cols: rho)")
    .print_matrix(res_random$mean_FDP, "mean_FDP (rows: kappa, cols: rho)")
  }

  # Run demo for Part 5
  trx_da_04_nn_05()

} # end Part 5
# ==============================================================================


# ==============================================================================
# Part 6: Monte Carlo simulation — 2D SNR x Rho sweep
# Fixed: n=300, p=1000, s=10, kappa=3, tFDR=0.2, K=20, 200 MC.
# Swept: rho in {0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9}
#        SNR in {0.1, 0.5, 1.0, 2.0, 5.0}
#
# Goal: identify whether higher signal strength can rescue the selector when
#       the correlation geometry is unfavourable.
#
# Compares CappedSpread (max_gap=20) vs Random support.
# ==============================================================================

if (TRUE) {

  trx_da_04_nn_06 <- function() {

    snr_grid <- c(0.1, 0.5, 1.0, 2.0, 5.0)
    rho_grid <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9)

    .run_snr_rho_sweep <- function(support_type, support_label) {

      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      n_s <- length(snr_grid)
      n_r <- length(rho_grid)

      row_nms <- paste0("SNR=", sprintf("%.2f", snr_grid))
      col_nms <- paste0("rho=", sprintf("%.2f", rho_grid))

      mat_TPP <- matrix(NA_real_, nrow = n_s, ncol = n_r,
                        dimnames = list(row_nms, col_nms))
      mat_FDP <- matrix(NA_real_, nrow = n_s, ncol = n_r,
                        dimnames = list(row_nms, col_nms))

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("NN MC — 2D SNR x Rho sweep\n")
      cat(sprintf("tFDR: %.2f | Support: %s | kappa=%d | %d MC\n",
                  MC_BASE$tFDR, support_label, MC_BASE$kappa, MC_BASE$num_MC))
      cat(" Sweeping SNR (rows) and rho (cols)\n")
      if (!is.null(support)) {
        cat(sprintf(" Support: {%s}\n", paste(support, collapse = ", ")))
      } else {
        cat(" Support: redrawn each trial\n")
      }
      cat(strrep("=", 70), "\n\n")
      # ------------------------------------------------

      total_cells <- n_s * n_r
      cell_idx <- 0L

      for (i_s in seq_len(n_s)) {
        snr_val <- snr_grid[i_s]

        for (i_r in seq_len(n_r)) {
          rho_val <- rho_grid[i_r]
          cell_idx <- cell_idx + 1L

          prefix <- sprintf(" [%2d/%2d] SNR=%-4.2f rho=%.2f",
                            cell_idx, total_cells, snr_val, rho_val)
          cat(sprintf("%s running %d MC trials ...\n", prefix, MC_BASE$num_MC))

          r <- .run_mc_nn(
            n = MC_BASE$n,
            p = MC_BASE$p,
            support = support,
            amplitude = MC_BASE$amplitude,
            kappa = MC_BASE$kappa,
            rho = rho_val,
            snr = snr_val,
            tFDR = MC_BASE$tFDR,
            K = MC_BASE$K,
            num_MC = MC_BASE$num_MC,
            seed = MC_BASE$seed,
            n_cores = num_cores
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
    res_capped <- .run_snr_rho_sweep(
      "capped", sprintf("CappedSpread(max_gap=%d)", MC_BASE$max_gap)
    )
    cat(sprintf("\n Results — CappedSpread(max_gap=%d):\n", MC_BASE$max_gap))
    .print_matrix(res_capped$mean_TPP, "mean_TPP (rows: SNR, cols: rho)")
    .print_matrix(res_capped$mean_FDP, "mean_FDP (rows: SNR, cols: rho)")

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

  # Run demo for Part 6
  trx_da_04_nn_06()

} # end Part 6
# ==============================================================================

cat("\nNN (MA(kappa)) demo complete.\n")
