# ==============================================================================
# demo_trex_spca_01_mc_sim.R
# ==============================================================================
#
# Monte Carlo simulation comparing T-Rex SPCA against PCA baselines.
# R counterpart of cpp demo_trex_spca_01_mc_sim.cpp.
#
# DGP: sparse M-factor model — X = Z * V^T + E.
#      n=50, p=100, M=3 latent factors, p1=5 active loadings per factor,
#      overlap_pool_size=30 (shared candidate-index pool across factors).
#
# Section 1: SNR sweep over the snr_values grid set below
#            (currently the single point -10 dB, num_MC=200).
#   Methods compared:
#     1. OrdPCA             — ordinary PCA (baseline; no sparsity constraint)
#     2. OraclePCA          — thresholded PCA with known support cardinality p1
#     3. TRexSPCA-EN-Act    — T-Rex+GVS(EN/TENET, Gram-based) + ActiveSet
#     4. TRexSPCA-ENaug-Act — T-Rex+GVS(EN/TENETAug, augmented) + ActiveSet
#     5. TRexSPCA-EN-Thr    — T-Rex+GVS(EN/TENET, Gram-based) + Thresholded
#     6. TRexSPCA-ENaug-Thr — T-Rex+GVS(EN/TENETAug, augmented) + Thresholded
#
#   All methods see the same center-only X (covariance-PCA footing, as in the
#   legacy R reference); the per-PC selector uses the default L2 scaling.
#
# Metrics (num_MC-averaged; TPR/FDR on PC1 support only — see evaluate_spca):
#   - FDR  — false discovery rate of PC1's estimated loading support
#   - TPR  — true positive rate  of PC1's estimated loading support
#   - PEV  — cumulative percentage of explained variance (all M components)
#
# Usage: Rscript demo_trex_spca_01_mc_sim.R [num_worker_cores]
# ==============================================================================

library(TRexSelectorNeo)
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

source(file.path(this_dir_, "..", "trex_spca_sim_utils.R"))


# ==============================================================================
# Monte Carlo Simulation — SNR sweep
# ==============================================================================

#' Run the SNR-sweep MC simulation and save results.
#'
#' @param cfg        List of simulation parameters (n, p, p1, M,
#'                   overlap_pool_size, tFDR, lambda_2, num_MC, base_seed).
#' @param snr_values SNR grid in dB.
#' @param stem       Output file stem (no folder, no extension).
demo_trex_spca_mc_snr_sweep <- function(cfg, snr_values, stem) {
  cat("================================================================\n")
  cat("  T-Rex SPCA MC Simulation — SNR Sweep\n")
  cat("================================================================\n")
  cat(sprintf("  n=%d  p=%d  M=%d  p1=%d  tFDR=%g  num_MC=%d  (%d cores)\n\n",
              cfg$n, cfg$p, cfg$M, cfg$p1, cfg$tFDR, cfg$num_MC, num_cores))

  method_names <- c("OrdPCA", "OraclePCA",
                    "TRexSPCA-EN-Act", "TRexSPCA-ENaug-Act",
                    "TRexSPCA-EN-Thr", "TRexSPCA-ENaug-Thr")

  empty_grid <- setNames(
    rep(list(numeric(length(snr_values))), length(method_names)), method_names)
  fdr_map <- tpr_map <- pev_map <- empty_grid

  for (si in seq_along(snr_values)) {
    store <- function(name, r) {
      fdr_map[[name]][si] <<- r$avg_fdr
      tpr_map[[name]][si] <<- r$avg_tpr
      pev_map[[name]][si] <<- r$avg_pev
    }
    snr <- snr_values[si]
    offset <- cfg$base_seed + as.integer(round(snr))
    snr_lbl <- sprintf("SNR=%.1f dB", snr)

    cat("--------------------------------------------------\n")
    cat(snr_lbl, "\n\n")

    make_data <- function(seed) {
      generate_sparse_factor_model(cfg$n, cfg$p, cfg$p1, cfg$M, snr,
                                   cfg$overlap_pool_size, seed)
    }

    # 1. Ordinary PCA
    r <- run_mc_trials_pca(cfg$num_MC, paste0(snr_lbl, " [OrdPCA]"),
                           make_data, offset, num_cores)
    store("OrdPCA", r)

    # 2. Oracle Thresholded PCA
    r <- run_mc_trials_oracle_pca(cfg$num_MC, paste0(snr_lbl, " [OraclePCA]"),
                                  make_data, cfg$p1, offset, num_cores)
    store("OraclePCA", r)

    # 3. T-Rex SPCA — EN (Gram-based TENET) + ActiveSet
    r <- run_mc_trials_trex_spca(cfg$num_MC, paste0(snr_lbl, " [TRexSPCA-EN-Act]"),
                                 make_data, cfg$tFDR, "ActiveSet", offset,
                                 cfg$lambda_2, "TENET", num_cores)
    store("TRexSPCA-EN-Act", r)

    # 4. T-Rex SPCA — ENaug (augmented TENET) + ActiveSet
    r <- run_mc_trials_trex_spca(cfg$num_MC, paste0(snr_lbl, " [TRexSPCA-ENaug-Act]"),
                                 make_data, cfg$tFDR, "ActiveSet", offset,
                                 cfg$lambda_2, "TENET_AUG", num_cores)
    store("TRexSPCA-ENaug-Act", r)

    # 5. T-Rex SPCA — EN (Gram-based TENET) + Thresholded
    r <- run_mc_trials_trex_spca(cfg$num_MC, paste0(snr_lbl, " [TRexSPCA-EN-Thr]"),
                                 make_data, cfg$tFDR, "Thresholded", offset,
                                 cfg$lambda_2, "TENET", num_cores)
    store("TRexSPCA-EN-Thr", r)

    # 6. T-Rex SPCA — ENaug (augmented TENET) + Thresholded
    r <- run_mc_trials_trex_spca(cfg$num_MC, paste0(snr_lbl, " [TRexSPCA-ENaug-Thr]"),
                                 make_data, cfg$tFDR, "Thresholded", offset,
                                 cfg$lambda_2, "TENET_AUG", num_cores)
    store("TRexSPCA-ENaug-Thr", r)
  }

  save_and_print_spca_mc_results(
    out_dir = file.path(this_dir_, "simulation_results"),
    file_stem = stem, num_MC = cfg$num_MC,
    snr_values = snr_values, method_names = method_names,
    fdr_map = fdr_map, tpr_map = tpr_map, pev_map = pev_map)
}


# ==============================================================================
# Main
# ==============================================================================

# =========================================================
# Section 1: SNR sweep
# =========================================================
if (TRUE) {
  cfg <- list(
    n                 = 50L,
    p                 = 100L,
    p1                = 5L,
    M                 = 3L,
    overlap_pool_size = 30L,
    tFDR              = 0.10,
    num_MC            = 200L,
    base_seed         = 42L,
    lambda_2          = -1.0  # -1 triggers k-fold CV for lambda_2 selection
  )

  snr_values <- c(-10.0)  # TEMP fast-validation

  demo_trex_spca_mc_snr_sweep(cfg, snr_values, "demo_trex_spca_01_mc_sim")
}
