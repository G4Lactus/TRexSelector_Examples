# ==============================================================================
# trex_spca_sim_utils.R
# ==============================================================================
#
# Shared simulation utilities for the T-Rex SPCA demo files. R counterpart of
# cpp/trex_selector_methods/trex_spca/trex_spca_sim_utils.hpp.
#
#   Data generation:
#     generate_sparse_factor_model()  — sparse M-factor model DGP
#                                       (returns list: X, Z, V, E, actual_snr_db)
#
#   Preprocessing:
#     center_columns()                — column centering (covariance-PCA footing)
#
#   Evaluation:
#     evaluate_spca()                 — PC1 TPR/FDR + cumulative PEV (QR-adjusted)
#
#   MC loop functions (parallel via parallel::mclapply):
#     run_mc_trials_trex_spca()       — MC loop for T-Rex SPCA
#     run_mc_trials_pca()             — MC loop for the ordinary PCA baseline
#     run_mc_trials_oracle_pca()      — MC loop for oracle thresholded PCA
#
#   Output:
#     save_and_print_spca_mc_results() — table (console + .txt) + tidy .csv
#
# Reproducibility pattern (mirrors cpp):
#   - Each trial's DATA is drawn from a deterministic seed
#     (base_seed_offset + mc*1000, disjoint 1000-wide bands per trial).
#   - The T-Rex selector runs with seed = -1 (hardware entropy per trial), so
#     the computer-generated dummies are independent across trials — required
#     for a valid Monte Carlo FDR estimate.
#
# Preprocessing note: X is CENTERED ONLY (covariance PCA, matching the legacy
# R reference and TRexSPCA's own internal convention). The column scales must
# NOT be normalized: the factor amplitude signal lives in the column variances,
# and z-scoring destroys it (measured at -10 dB: T-Rex SPCA FDR degrades from
# ~0.14 to ~0.5 — see cpp validation_trex_spca_06_handrolled_comparison).
# ==============================================================================

library(TRexSelectorNeo)
library(parallel)


# ==============================================================================
# Data generator (sparse M-factor model)
# ==============================================================================

#' Generate a single dataset from the sparse factor model X = Z V^T + E.
#'
#' Z (n x M) has column-wise N(0, sigma_m^2) factors with sigma = {5, 3, 1}
#' (falling back to 1 for M > 3). V (p x M) has exactly p1 non-zeros (value 0.9)
#' per column, drawn without replacement from a shared candidate pool of size
#' overlap_pool_size — so different factors' active sets can overlap. E is
#' i.i.d. Gaussian noise scaled to match the target SNR in dB.
#'
#' @param n                  Number of samples.
#' @param p                  Number of features.
#' @param p1                 True active loadings per factor.
#' @param M                  Number of latent factors.
#' @param target_snr         Target signal-to-noise ratio in dB.
#' @param overlap_pool_size  Shared candidate-index pool size.
#' @param seed               RNG seed for this dataset.
#'
#' @return list(X, Z, V, E, actual_snr_db); X is raw (uncentered).
generate_sparse_factor_model <- function(n = 50L, p = 100L, p1 = 5L, M = 3L,
                                         target_snr = 0.0,
                                         overlap_pool_size = 30L,
                                         seed = 42L) {
  stopifnot(p1 <= p, n > 0L, p > 0L, p1 > 0L)
  set.seed(seed)

  # 1. True factors Z (n x M) with decreasing signal strength
  factor_stds <- c(5.0, 3.0, 1.0)
  Z <- matrix(0, nrow = n, ncol = M)
  for (m in seq_len(M)) {
    std_dev <- if (m <= length(factor_stds)) factor_stds[m] else 1.0
    Z[, m] <- rnorm(n, mean = 0, sd = std_dev)
  }

  # 2. Sparse loadings V (p x M): each factor draws p1 indices from the pool
  actual_pool <- min(overlap_pool_size, p)
  stopifnot(p1 <= actual_pool)
  V <- matrix(0, nrow = p, ncol = M)
  for (m in seq_len(M)) {
    active_idx <- sample.int(actual_pool, p1, replace = FALSE)
    V[active_idx, m] <- 0.9
  }

  # 3. Signal matrix S = Z V^T
  S <- Z %*% t(V)
  signal_var <- stats::var(as.vector(S))

  # 4. Gaussian noise E scaled to the target SNR
  # SNR_dB = 10 * log10(var_S / var_E)  =>  var_E = var_S / 10^(SNR_dB / 10)
  E_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
  raw_noise_var <- stats::var(as.vector(E_raw))
  target_noise_var <- signal_var / (10^(target_snr / 10))
  E <- E_raw * sqrt(target_noise_var / raw_noise_var)

  # 5. Assemble X
  X <- S + E
  final_noise_var <- stats::var(as.vector(E))
  actual_snr_db <- 10 * log10(signal_var / final_noise_var)

  list(X = X, Z = Z, V = V, E = E, actual_snr_db = actual_snr_db)
}


# ==============================================================================
# Preprocessing
# ==============================================================================

#' Center columns (mean subtraction only, NO scaling).
#'
#' Covariance-PCA footing shared by every method, matching the legacy R
#' reference pipeline and TRexSPCA's internal convention. Column scales carry
#' the factor amplitude signal and must not be normalized.
center_columns <- function(X) {
  sweep(X, 2L, colMeans(X), "-")
}


# ==============================================================================
# Single-trial evaluator
# ==============================================================================

#' Compute per-trial sparse PCA metrics from estimated loadings and scores.
#'
#' TPR / FDR are evaluated strictly on the FIRST principal component (PC1):
#' PCA's orthogonality constraint mixes the loading supports of PCs 2+ across
#' the true factors, making per-component support recovery metrics ambiguous
#' beyond the first PC. Cumulative PEV is computed via QR decomposition of
#' Z_est across all M components.
#'
#' @param X       Centered observation matrix (same X the methods saw).
#' @param V_est   Estimated sparse loading matrix (p x M).
#' @param Z_est   Estimated score matrix (n x M).
#' @param V_true  True loading matrix (p x M); non-zeros mark the true support.
#'
#' @return list(tpr, fdr, pev)
evaluate_spca <- function(X, V_est, Z_est, V_true) {
  n <- nrow(X)
  M <- ncol(V_est)

  # 1. TPR / FDR strictly on PC1
  true_supp_pc1 <- which(V_true[, 1] != 0)
  est_supp_pc1  <- which(abs(V_est[, 1]) > 1e-12)

  tp <- length(intersect(true_supp_pc1, est_supp_pc1))
  tpr <- if (length(true_supp_pc1) == 0L) 1.0 else tp / length(true_supp_pc1)
  fdr <- if (length(est_supp_pc1) == 0L) 0.0 else
    (length(est_supp_pc1) - tp) / length(est_supp_pc1)

  # 2. Cumulative PEV across ALL components (QR of Z_est)
  R <- qr.R(qr(Z_est))
  total_var <- sum(X^2) / (n - 1)
  cum_ev <- sum(diag(R)[seq_len(M)]^2) / (n - 1)
  pev <- if (total_var > 0) cum_ev / total_var else 0.0

  list(tpr = tpr, fdr = fdr, pev = pev)
}


# ==============================================================================
# MC aggregation helper
# ==============================================================================

.aggregate_spca_mc <- function(label, tpr_vec, fdr_vec, pev_vec) {
  res <- list(
    avg_tpr = mean(tpr_vec), avg_fdr = mean(fdr_vec), avg_pev = mean(pev_vec),
    sd_tpr = if (length(tpr_vec) > 1L) stats::sd(tpr_vec) else 0.0,
    sd_fdr = if (length(fdr_vec) > 1L) stats::sd(fdr_vec) else 0.0
  )
  cat(sprintf("  %s — done. TPR=%.3f  FDR=%.3f\n\n", label,
              res$avg_tpr, res$avg_fdr))
  res
}

.mc_map <- function(num_MC, one_trial, num_cores) {
  # mc index runs 0..num_MC-1 (matches the cpp seed bands)
  res_list <- parallel::mclapply(0:(num_MC - 1L), one_trial,
                                 mc.cores = num_cores, mc.preschedule = FALSE)
  errs <- vapply(res_list, inherits, logical(1), what = "try-error")
  if (any(errs)) stop("MC trial failed: ", res_list[[which(errs)[1]]])
  res_list
}


# ==============================================================================
# Parallel MC inner loops
# ==============================================================================

#' Run num_MC parallel T-Rex SPCA trials.
#'
#' Per trial: make_data(seed) returns raw data; the pipeline centers X once
#' (covariance PCA); this same X is shared by the selector and by
#' evaluate_spca(). The selector runs with seed = -1 (hardware entropy) so
#' dummies are independent across trials.
#'
#' @param num_MC            Number of Monte Carlo trials.
#' @param progress_label    Label printed in start/done lines.
#' @param make_data         Factory: seed -> data list (called per trial).
#' @param tFDR              Target FDR for T-Rex SPCA.
#' @param mode              "ActiveSet" or "Thresholded".
#' @param base_seed_offset  Base seed; trial mc uses base_seed_offset + mc*1000.
#' @param lambda_2          Ridge penalty (LARS units); < 0 = auto via CV,
#'                          0 = no ridge, > 0 = fixed.
#' @param en_solver         "TENET" or "TENET_AUG".
#' @param num_cores         Worker processes for mclapply.
#'
#' @return list(avg_tpr, avg_fdr, avg_pev, sd_tpr, sd_fdr)
run_mc_trials_trex_spca <- function(num_MC, progress_label, make_data,
                                    tFDR, mode, base_seed_offset,
                                    lambda_2 = -1.0, en_solver = "TENET",
                                    num_cores = 1L) {
  cat(sprintf("  %s — Running %d MC trials ...\n", progress_label, num_MC))

  one_trial <- function(mc) {
    dat <- make_data(base_seed_offset + mc * 1000L)
    X <- center_columns(dat$X)
    M <- ncol(dat$V)

    # X + 0 forces a fresh copy: the selector maps X zero-copy via Eigen and
    # centers it internally; the local X must stay intact for evaluation.
    sel <- TRexSPCASelector$new(
      X + 0, M = M, tFDR = tFDR,
      control = trex_spca_control(mode = mode, gvs_type = "EN",
                                  en_solver = en_solver,
                                  lambda_2 = lambda_2,
                                  lambda2_method = "CV_1SE_CCD"),
      seed = -1L, verbose = FALSE)
    sel$select()

    unlist(evaluate_spca(X, sel$V, sel$Z, dat$V))
  }

  res_mat <- do.call(rbind, .mc_map(num_MC, one_trial, num_cores))
  .aggregate_spca_mc(progress_label,
                     res_mat[, "tpr"], res_mat[, "fdr"], res_mat[, "pev"])
}


#' Run num_MC parallel ordinary PCA trials.
#'
#' Ordinary PCA retains all p loadings per component (no sparsity), so
#' FDR ~ (p - p1) / p and TPR = 1 regardless of SNR; the metric of interest
#' is the cumulative PEV.
run_mc_trials_pca <- function(num_MC, progress_label, make_data,
                              base_seed_offset, num_cores = 1L) {
  cat(sprintf("  %s — Running %d MC trials ...\n", progress_label, num_MC))

  one_trial <- function(mc) {
    dat <- make_data(base_seed_offset + mc * 1000L)
    X <- center_columns(dat$X)
    M <- ncol(dat$V)

    V_ord <- svd(X, nu = 0, nv = M)$v
    Z_ord <- X %*% V_ord

    unlist(evaluate_spca(X, V_ord, Z_ord, dat$V))
  }

  res_mat <- do.call(rbind, .mc_map(num_MC, one_trial, num_cores))
  .aggregate_spca_mc(progress_label,
                     res_mat[, "tpr"], res_mat[, "fdr"], res_mat[, "pev"])
}


#' Run num_MC parallel oracle-thresholded PCA trials.
#'
#' Oracle PCA keeps exactly the top p1 entries (by absolute value) of each
#' ordinary-PCA loading vector — it knows the true support cardinality. This
#' bounds what any data-driven thresholding method can achieve.
run_mc_trials_oracle_pca <- function(num_MC, progress_label, make_data,
                                     p1, base_seed_offset, num_cores = 1L) {
  cat(sprintf("  %s — Running %d MC trials ...\n", progress_label, num_MC))

  one_trial <- function(mc) {
    dat <- make_data(base_seed_offset + mc * 1000L)
    X <- center_columns(dat$X)
    p <- ncol(X)
    M <- ncol(dat$V)

    V_ord <- svd(X, nu = 0, nv = M)$v

    # Keep top p1 entries per loading vector, then normalize
    V_oracle <- matrix(0, nrow = p, ncol = M)
    for (comp in seq_len(M)) {
      v_m <- V_ord[, comp]
      thresh <- sort(abs(v_m), decreasing = TRUE)[p1]
      keep <- abs(v_m) >= thresh
      V_oracle[keep, comp] <- v_m[keep]
      V_oracle[, comp] <- V_oracle[, comp] / sqrt(sum(V_oracle[, comp]^2))
    }
    Z_oracle <- X %*% V_oracle

    unlist(evaluate_spca(X, V_oracle, Z_oracle, dat$V))
  }

  res_mat <- do.call(rbind, .mc_map(num_MC, one_trial, num_cores))
  .aggregate_spca_mc(progress_label,
                     res_mat[, "tpr"], res_mat[, "fdr"], res_mat[, "pev"])
}


# ==============================================================================
# Output
# ==============================================================================

#' Print MC-averaged T-Rex SPCA results as an aligned table (console + .txt)
#' and write a tidy long-format CSV (method, metric, snr_db, value).
#'
#' @param out_dir       Output directory (created if missing).
#' @param file_stem     Output file base name (no folder, no extension).
#' @param num_MC        Number of Monte Carlo trials.
#' @param snr_values    SNR grid (dB) used during the simulation.
#' @param method_names  Ordered vector of method names (keys into the maps).
#' @param fdr_map       Named list: method -> numeric vector over the SNR grid.
#' @param tpr_map       Named list: method -> numeric vector over the SNR grid.
#' @param pev_map       Named list: method -> numeric vector over the SNR grid.
save_and_print_spca_mc_results <- function(out_dir, file_stem, num_MC,
                                           snr_values, method_names,
                                           fdr_map, tpr_map, pev_map) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  txt_path <- file.path(out_dir, paste0(file_stem, ".txt"))
  con <- file(txt_path, open = "wt")
  sink(con, split = TRUE)
  on.exit({ sink(); close(con) }, add = TRUE)

  cat("\n")
  cat("======================================================================\n")
  cat(sprintf("=== T-Rex SPCA Results (averaged over %d MC trials) ===\n", num_MC))
  cat("======================================================================\n\n")

  # Column header + dashed separator
  hdr <- sprintf("%-22s%-8s%8s", "Method", "Metric", "SNR(dB)")
  hdr <- paste0(hdr, paste(sprintf("%10.1f", snr_values), collapse = ""))
  cat(hdr, "\n", sep = "")
  cat(strrep("-", 38L + 10L * length(snr_values)), "\n", sep = "")

  print_row <- function(mname, metric, values, first_row) {
    cat(sprintf("%-22s%-8s%8s", if (first_row) mname else "", metric, ""))
    cat(sprintf("%10.4f", values), sep = "")
    cat("\n")
  }

  for (name in method_names) {
    print_row(name, "FDR", fdr_map[[name]], TRUE)
    print_row(name, "TPR", tpr_map[[name]], FALSE)
    print_row(name, "PEV", pev_map[[name]], FALSE)
    cat("\n")
  }

  sink()
  close(con)
  on.exit(NULL)

  # Tidy long-format CSV (method, metric, snr_db, value)
  rows <- list()
  for (name in method_names) {
    for (i in seq_along(snr_values)) {
      rows[[length(rows) + 1L]] <- data.frame(
        method = name,
        metric = c("FDR", "TPR", "PEV"),
        snr_db = snr_values[i],
        value  = c(fdr_map[[name]][i], tpr_map[[name]][i], pev_map[[name]][i]),
        stringsAsFactors = FALSE)
    }
  }
  csv_path <- file.path(out_dir, paste0(file_stem, ".csv"))
  utils::write.csv(do.call(rbind, rows), csv_path,
                   row.names = FALSE, quote = FALSE)

  cat(sprintf("[Info] CSV results saved to:             %s\n", csv_path))
  cat(sprintf("[Info] Results successfully saved to: %s\n\n", txt_path))
}
