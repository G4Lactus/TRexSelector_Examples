# ==============================================================================
# demo_trex_da_01_ar1.R
# ==============================================================================
#
# DA-TRex demo and MC simulation for the AR(1) Toeplitz DGP.
#
#  Part 1: Single-run demo.
#
#  Part 2: MC simulation sweeping SNR. Compares Random support vs.
#          CappedSpread (max_gap=20).
#
#  Part 3: MC simulation sweeping rho. CappedSpread support (max_gap=20).
#
#  Part 4: 2D Gap x Rho sweep. Explores the kappa-boundary where gap < kappa
#          collapses TPP.
#
# ==============================================================================

# Load packages
library(TRexSelector)
library(plotly)
library(parallel)

num_cores <- 6

# Resolve the directory of this file, works for both source() and Rscript
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) "."
)

# Source files
source(file.path(this_dir_, "..", "support_generators.R"))
source(file.path(this_dir_, "dgp_ar1_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))


# ==============================================================================
# Base parameters
# ==============================================================================

# Part 1 relevant parameters
PARAMS <- list(
  n         = 150L,   # observations
  p         = 500L,   # predictors
  s         = 5L,     # active predictors
  amplitude = 3.0,    # signal coefficient for each active beta_j
  snr       = 2.0,    # signal-to-noise ratio
  K         = 20L,    # random experiments per selector run
  tFDR      = 0.2,    # target FDR level
  seed      = 2026L,  # global RNG seed
  rho       = 0.7,    # AR(1) correlation
  max_gap   = 15L     # max gap for CappedSpread support
)

# Parts 2-4: Monte Carlo simulation base parameters
MC_BASE <- list(
  n         = 300L,
  p         = 1000L,
  s         = 10L,
  amplitude = 3.0,
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
  cat(sprintf(" Data: n = %d, p = %d, tFDR = %.2f, s = %d\n",
              dat$n, dat$p, tFDR, dat$s))
  cat(sprintf("       beta_active = %.2f, SNR = %.2f,  rho = %.2f\n",
              PARAMS$amplitude, dat$snr, PARAMS$rho))
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
# Part 1: Single-run demo
# ==============================================================================

if (FALSE) {
  trx_da_01_ar1_01 <- function() {

    # Header
    # -------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("AR(1) Demo  |  Part: 1 Single-run\n")
    cat(sprintf("n=%d,  p=%d,  s=%d,  rho=%.2f, tFDR=%.2f, beta_active=%.2f\n",
      PARAMS$n, PARAMS$p, PARAMS$s, PARAMS$rho, PARAMS$tFDR, PARAMS$amplitude
    ))
    cat(strrep("=", 70), "\n\n")
    # -------------------------------------------------

    # Support set generation
    random_indexing <- TRUE
    if (random_indexing) {
      # Random support sampling
      support_p1 <- make_support_random(PARAMS$s, PARAMS$p, PARAMS$seed)
    } else {
      # CappedSpread: deterministic, evenly-spaced
      support_p1 <- make_support_capped_spread(PARAMS$s, PARAMS$p, PARAMS$max_gap)
    }

    # Generate data with SNR-calibrated noise
    cat("[Part 1] Generating AR(1) data ...\n")
    dat_p1 <- dgp_ar1_snr(
      n         = PARAMS$n,
      p         = PARAMS$p,
      support   = support_p1,
      amplitude = PARAMS$amplitude,
      rho       = PARAMS$rho,
      snr       = PARAMS$snr,
      seed      = PARAMS$seed
    )

    # Visualize the correlation matrix of the predictors
    cor_mat_d_AR1 <- cor(dat_p1$X)
    plot_cormat(cor_mat_d_AR1)

    # Run T-Rex+DA+AR1 (LARS solver, auto-estimated rho)
    cat("[Part 1] Running trex+DA+AR1 (cor_coef = NA -> auto-estimated) ...\n\n")
    res_p1 <- TRexSelector::trex(
      X          = dat_p1$X,
      y          = dat_p1$y,
      tFDR       = PARAMS$tFDR,
      K          = PARAMS$K,
      method     = "trex+DA+AR1",
      cor_coef   = NA,
      type       = "lar",
      rho_thr_DA = 0.02,
      verbose    = FALSE,
      seed       = PARAMS$seed
    )

    .print_result(
      "Part 1 — AR(1) Demo  [trex+DA+AR1, TLARS, cor_coef=auto]",
      dat_p1, res_p1, PARAMS$tFDR
    )
  }

  # Run demo for Part 1
  trx_da_01_ar1_01()

}  # end Part 1
# ==============================================================================



# ==============================================================================
# Part 2: Monte Carlo simulation — SNR sweep
# ==============================================================================
# Compares Random support (redrawn per trial) vs CappedSpread support.
# Fixed: n=300, p=1000, s=10, amplitude=3.0, rho=0.7, tFDR=0.2, K=20, 200 MC.
# SNR grid: {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}.
# ==============================================================================

if (FALSE) {
  trx_da_01_ar1_02 <- function() {

    snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

    # Helper: run SNR sweep for one support type.
    .run_snr_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL  # redrawn per trial inside .run_mc_ar1()
      }

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("AR(1) MC — SNR sweep\n")
      cat(sprintf("tFDR: %.2f  |  Support: %s  |  amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d  p=%d  s=%d  rho=%.1f  SNR swept  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$rho, MC_BASE$num_MC))
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
        r <- .run_mc_ar1(MC_BASE$n, MC_BASE$p, support, MC_BASE$amplitude,
                         MC_BASE$rho, snr_val,
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
  trx_da_01_ar1_02()

}  # end Part 2
# ==============================================================================



# ==============================================================================
#  Part 3: Monte Carlo simulation — Rho sweep
#
#  Fixed: n=300, p=1000, s=10, amplitude=3.0, SNR=2.0, tFDR=0.2, K=20, 200 MC.
#  Compares CappedSpread (max_gap=20) vs Random support; rho in {0.0, ..., 0.9}.
# ==============================================================================

if (FALSE) {
  trx_da_01_ar1_03 <- function() {

    rho_grid <- seq(0.0, 0.9, by = 0.1)

    # Run rho sweep for one support type.
    .run_rho_sweep <- function(support_type, support_label) {
      support <- if (support_type == "capped") {
        make_support_capped_spread(MC_BASE$s, MC_BASE$p, MC_BASE$max_gap)
      } else {
        NULL
      }

      # Header
      # ------------------------------------------------
      cat(strrep("=", 70), "\n")
      cat("AR(1) MC — Rho sweep\n")
      cat(sprintf("tFDR: %.2f  |  Support: %s  |  amplitude: %.2f\n",
                  MC_BASE$tFDR, support_label, MC_BASE$amplitude))
      cat(sprintf("n=%d  p=%d  s=%d  rho swept  SNR=%.1f  %d MC\n",
                  MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$snr, MC_BASE$num_MC))
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
        r <- .run_mc_ar1(MC_BASE$n, MC_BASE$p, support, MC_BASE$amplitude,
                         rho_val, MC_BASE$snr,
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
  trx_da_01_ar1_03()

}  # end Part 3
# ==============================================================================


# ==============================================================================
#  Part 4: 2D sweep: Gap x Rho
#
#  Runs 200 MC trials per (gap, rho) cell for two support types:
#    A) CappedSpread(gap) — deterministic, evenly-spaced support.
#    B) Random           — support redrawn each trial (not gap-dependent).
#
#  The kappa boundary (kappa = ceiling(log(0.02) / log(rho))) is the DA+AR1
#  correction window half-width.  When gap < kappa, active predictors fall
#  inside each other's correction windows and TPP collapses.
#
#  gap_grid : {100, 50, 20, 15, 10, 5, 1}
#  rho_grid : {0.0, 0.1, ..., 0.9}
# ==============================================================================

if (FALSE) {
  trx_da_01_ar1_04 <- function() {

    gap_grid <- c(100L, 50L, 20L, 15L, 10L, 5L, 1L)
    rho_grid <- seq(0.0, 0.9, by = 0.1)

    # kappa: smallest window half-width such that rho^kappa <= 0.02
    kappa_boundary <- ceiling(log(0.02) / log(rho_grid))

    n_gap <- length(gap_grid)
    n_rho <- length(rho_grid)

    row_nms <- paste0("rho=", sprintf("%.1f", rho_grid))
    col_nms <- paste0("gap=", gap_grid)

    mat_TPP_cs <- matrix(NA_real_, nrow = n_rho, ncol = n_gap,
                         dimnames = list(row_nms, col_nms))
    mat_FDP_cs <- mat_TPP_cs
    rand_TPP   <- setNames(numeric(n_rho), row_nms)
    rand_FDP   <- setNames(numeric(n_rho), row_nms)

    # Header
    # ------------------------------------------------
    cat(strrep("=", 70), "\n")
    cat("AR(1) Demo  |  Part: 4 Gap x Rho sweep\n")
    cat(sprintf(
      "n=%d, p=%d, s=%d, SNR=%.1f, %d MC, tFDR=%.2f, amplitude=%.2f\n",
      MC_BASE$n, MC_BASE$p, MC_BASE$s, MC_BASE$snr, MC_BASE$num_MC,
      MC_BASE$tFDR, MC_BASE$amplitude
    ))
    cat(sprintf("  gap_grid  : {%s}\n", paste(gap_grid, collapse = ", ")))
    cat(sprintf("  rho_grid  : {%s}\n", paste(sprintf("%.1f", rho_grid), collapse = ", ")))
    cat(strrep("=", 70), "\n\n")
    # ------------------------------------------------


    # ---- Section A: CappedSpread sweep (gap x rho) ----
    cat("  [A] CappedSpread(gap) sweep ...\n\n")

    total_cells <- n_rho * n_gap
    cell_idx    <- 0L

    for (i_r in seq_len(n_rho)) {
      rho_val <- rho_grid[i_r]
      for (i_g in seq_len(n_gap)) {
        gap_val  <- gap_grid[i_g]
        cell_idx <- cell_idx + 1L

        support_mc <- make_support_capped_spread(
          s       = MC_BASE$s,
          p       = MC_BASE$p,
          max_gap = gap_val
        )

        prefix <- sprintf(
          "  [%d/%d] rho=%.1f gap=%3d kappa=%-3d",
          cell_idx, total_cells, rho_val, gap_val,
          kappa_boundary[i_r]
        )
        cat(sprintf("%s running %d MC trials ...\n", prefix, MC_BASE$num_MC))

        r <- .run_mc_ar1(MC_BASE$n, MC_BASE$p, support_mc, MC_BASE$amplitude,
                         rho_val, MC_BASE$snr,
                         MC_BASE$tFDR, MC_BASE$K, MC_BASE$num_MC, MC_BASE$seed,
                         n_cores = num_cores)

        mat_TPP_cs[i_r, i_g] <- r$mean_TPP
        mat_FDP_cs[i_r, i_g] <- r$mean_FDP
        cat(sprintf("    %s done. TPP=%.3f FDP=%.3f\n\n",
                    prefix, r$mean_TPP, r$mean_FDP))
      }
    }


    # ---- Section B: Random support sweep (rho only) ----
    cat("\n  [B] Random support sweep (support redrawn each trial) ...\n\n")

    for (i_r in seq_len(n_rho)) {
      rho_val <- rho_grid[i_r]

      prefix <- sprintf(
        "  [%d/%d] rho=%.1f kappa=%-3d Random support",
        i_r, n_rho, rho_val, kappa_boundary[i_r]
      )
      cat(sprintf("%s running %d MC trials ...\n", prefix, MC_BASE$num_MC))

      r <- .run_mc_ar1(MC_BASE$n, MC_BASE$p, NULL, MC_BASE$amplitude,
                       rho_val, MC_BASE$snr,
                       MC_BASE$tFDR, MC_BASE$K, MC_BASE$num_MC, MC_BASE$seed,
                       n_cores = num_cores)

      rand_TPP[i_r] <- r$mean_TPP
      rand_FDP[i_r] <- r$mean_FDP
      cat(sprintf("    %s done. TPP=%.3f FDP=%.3f\n\n",
                  prefix, r$mean_TPP, r$mean_FDP))
    }


    # ---- Print results ----

    # kappa reference table
    cat("\n  DA+AR1 correction window half-width: kappa = ceiling(log(0.02) / log(rho))\n")
    cat(sprintf("  %-8s  %s\n", "rho",
                paste(sprintf("%-5.1f", rho_grid), collapse = "  ")))
    cat(sprintf("  %-8s  %s\n", "kappa",
                paste(sprintf("%-5d", kappa_boundary), collapse = "  ")))

    # Combine CappedSpread matrix with Random column for side-by-side comparison
    mat_TPP_combined <- cbind(mat_TPP_cs, Random = rand_TPP)
    mat_FDP_combined <- cbind(mat_FDP_cs, Random = rand_FDP)

    .print_matrix(mat_TPP_combined, "mean_TPP  (cols: CappedSpread gap values + Random)")
    .print_matrix(mat_FDP_combined, "mean_FDP  (cols: CappedSpread gap values + Random)")
  }

  # Run demo for Part 4
  trx_da_01_ar1_04()

}  # end Part 4
# ==============================================================================

cat("\nAR(1) demo complete.\n")
