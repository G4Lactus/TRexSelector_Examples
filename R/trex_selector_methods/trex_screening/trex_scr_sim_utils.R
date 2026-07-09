# ==============================================================================
# trex_scr_sim_utils.R
# ==============================================================================
#
# Shared infrastructure for the Screen-TRex R demonstration suite.
#
# Mirrors the two C++ demo-internal headers:
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_common.hpp
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_bb_common.hpp
#
# Provides:
#   Section 1 — data-generating processes (i.i.d., AR(1), equi, block-equi),
#   Section 2 — Screen-TRex control-parameter factories,
#   Section 3 — single-run result printer,
#   Section 4 — Monte Carlo runner + dual TXT/CSV table output (Screen-TRex),
#   Section 5 — Biobank (Algorithm 1) helpers + MC table output.
#
# NOTE on reproducibility: the C++ demos draw with std::mt19937; R uses its own
# RNG stream. Seeds are mirrored as *labels*, but the realised random draws
# differ, so R results are statistically comparable to the C++ reference — not
# bit-identical.
#
# Support sets and selected indices are 1-BASED throughout (R convention); the
# C++ sources are 0-based.
# ==============================================================================

library(TRexSelectorNeo)


# ==============================================================================
# Section 1 — Data-generating processes
# ==============================================================================

#' Evenly-spaced sparse beta vector (mirrors make_beta() in demo 03).
#'
#' Active indices are placed at round((i)*p/(p1+2)) for i = 1..p1 (0-based in
#' C++), converted to 1-based R indices.
#'
#' @param p        Total number of predictors.
#' @param p1       Number of active predictors.
#' @param beta_val Coefficient value for active predictors.
#' @return Numeric vector of length p with p1 non-zero, evenly spaced entries.
make_beta <- function(p, p1, beta_val) {
  beta <- numeric(p)
  for (i in seq_len(p1)) {
    idx0 <- round(i * p / (p1 + 2))       # 0-based position (C++)
    beta[idx0 + 1L] <- beta_val            # 1-based (R)
  }
  beta
}

#' SNR-calibrated response y = X beta + eps, then mean-centred.
#' Mirrors make_response(): sigma = sqrt(||X beta||^2 / (n * snr)).
.make_response <- function(X, beta, snr) {
  n      <- nrow(X)
  signal <- as.numeric(X %*% beta)
  sigma  <- sqrt(sum(signal^2) / (n * snr))
  y      <- signal + rnorm(n, 0, sigma)
  y - mean(y)
}

#' Standardise columns to mean 0 and unit variance (n - 1 denominator).
.standardise_cols <- function(X) {
  X <- sweep(X, 2L, colMeans(X), "-")
  sds <- sqrt(colSums(X^2) / (nrow(X) - 1L))
  sds[sds < 1e-12] <- 1
  sweep(X, 2L, sds, "/")
}

#' i.i.d. Gaussian DGP (mirrors datagen::SyntheticData used in demos 01/02/06).
#'
#' X ~ N(0, 1) i.i.d.; y = X beta + eps with SNR-calibrated noise.
#'
#' @param n,p          Dimensions.
#' @param true_support 1-based indices of active variables.
#' @param coefs        Coefficients for the active variables.
#' @param snr          Signal-to-noise ratio.
#' @param seed         RNG seed.
#' @return list(X, y, true_support).
dgp_iid_snr <- function(n, p, true_support, coefs, snr, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[true_support] <- coefs
  y <- .make_response(X, beta, snr)
  list(X = X, y = y, true_support = sort(true_support))
}

#' AR(1) column-correlated DGP (mirrors dgp_ar1 in demo 03).
dgp_ar1 <- function(n, p, p1, rho, snr, beta_val, seed) {
  set.seed(seed)
  innov_sd <- sqrt(1 - rho^2)
  X <- matrix(0, n, p)
  X[, 1L] <- rnorm(n)
  for (j in 2:p) X[, j] <- rho * X[, j - 1L] + innov_sd * rnorm(n)
  X <- .standardise_cols(X)
  beta <- make_beta(p, p1, beta_val)
  y <- .make_response(X, beta, snr)
  list(X = X, y = y, true_support = which(beta != 0))
}

#' Equicorrelated DGP (mirrors dgp_equi in demo 03).
#' x_j = sqrt(rho) * f + sqrt(1 - rho) * w_j, shared factor f per observation.
dgp_equi <- function(n, p, p1, rho, snr, beta_val, seed) {
  set.seed(seed)
  load_f   <- sqrt(rho)
  load_eta <- sqrt(1 - rho)
  f <- rnorm(n)
  X <- load_f * matrix(f, n, p) + load_eta * matrix(rnorm(n * p), n, p)
  X <- .standardise_cols(X)
  beta <- make_beta(p, p1, beta_val)
  y <- .make_response(X, beta, snr)
  list(X = X, y = y, true_support = which(beta != 0))
}

#' Block-equicorrelated DGP (mirrors dgp_block_equi in demo 03).
#' Contiguous blocks, each with its own shared factor.
dgp_block_equi <- function(n, p, p1, rho, n_blocks, snr, beta_val, seed) {
  set.seed(seed)
  load_f    <- sqrt(rho)
  load_eta  <- sqrt(1 - rho)
  base_size <- p %/% n_blocks
  remainder <- p %% n_blocks
  X <- matrix(0, n, p)
  col <- 0L
  for (k in seq_len(n_blocks)) {
    bk <- base_size + if (k <= remainder) 1L else 0L
    z_k <- rnorm(n)                                   # shared factor for block k
    idx <- (col + 1L):(col + bk)
    X[, idx] <- load_f * matrix(z_k, n, bk) +
                load_eta * matrix(rnorm(n * bk), n, bk)
    col <- col + bk
  }
  X <- .standardise_cols(X)
  beta <- make_beta(p, p1, beta_val)
  y <- .make_response(X, beta, snr)
  list(X = X, y = y, true_support = which(beta != 0))
}


# ==============================================================================
# Section 2 — Control-parameter factories
# ==============================================================================

#' Default TRexControlParameter for Screen-TRex (mirrors make_trex_control()).
#'
#' @param solver             T-Rex solver backend (default "TLARS").
#' @param use_memory_mapping Enable D-mmap + solver serialization.
#' @param ...                Extra trex_control() overrides (e.g. rho_afs).
make_scr_trex_ctrl <- function(solver = "TLARS",
                               use_memory_mapping = FALSE, ...) {
  trex_control(
    solver               = solver,
    K                    = 20L,
    max_dummy_multiplier = 10L,
    use_max_T_stop       = TRUE,
    lloop_strategy       = "STANDARD",
    use_memory_mapping   = use_memory_mapping,
    ...
  )
}

#' Screen-TRex method descriptors for the Ordinary-vs-Bootstrap comparison
#' (mirrors default_methods()).
default_scr_methods <- function() {
  list(
    list(name = "Screen-TRex Ordinary",  trex_method = "TREX", bootstrap = FALSE),
    list(name = "Screen-TRex Bootstrap", trex_method = "TREX", bootstrap = TRUE)
  )
}


# ==============================================================================
# Section 3 — Single-run result printer
# ==============================================================================

#' Print a single Screen-TRex selection block (mirrors print_screen_result()).
#'
#' @param sel          A selected TRexScreeningSelector object.
#' @param true_support 1-based ground-truth active indices.
.print_scr_result <- function(sel, true_support) {
  selected <- sel$selected_indices
  cat(sprintf("\n  Estimated FDR:  %.4f\n", sel$estimated_FDR))
  if (!is.null(sel$estimated_correlation) &&
      abs(sel$estimated_correlation) > 0) {
    cat(sprintf("  Estimated rho:  %.4f\n", sel$estimated_correlation))
  }
  if (isTRUE(sel$used_bootstrap)) {
    cat("  Bootstrap:      yes\n")
    cat(sprintf("  Conf. level:    %.4f\n", sel$confidence_level))
  }
  cat(sprintf("  Selected (%d):  %s\n",
              length(selected), paste(sort(selected), collapse = " ")))
  fdp <- compute_fdp(selected, true_support)
  tpp <- compute_tpp(selected, true_support)
  cat(sprintf("  FDP: %.4f  |  TPP: %.4f\n", fdp, tpp))
  cat(sprintf("  True support:   %s\n\n", paste(sort(true_support), collapse = " ")))
  invisible(list(fdp = fdp, tpp = tpp))
}


# ==============================================================================
# Section 4 — Monte Carlo runner + table output (Screen-TRex)
# ==============================================================================

#' Run num_MC Screen-TRex trials for one grid point.
#'
#' Each trial uses the DGP factory with seed = base_seed + mc - 1 and a
#' selector seed of -1 (hardware entropy per trial, mirroring the C++ MC loop).
#'
#' @param num_MC       Monte Carlo trials.
#' @param label        Progress label.
#' @param make_data    function(seed) -> list(X, y, true_support).
#' @param trex_ctrl    trex_control() object (constant across trials).
#' @param screen_args  named list of trex_screen_control() arguments.
#' @param base_seed    Trial seed base.
#' @return list(avg_fdr, avg_tpr, avg_est_fdr).
run_mc_screen <- function(num_MC, label, make_data, trex_ctrl,
                          screen_args, base_seed) {
  fdp <- tpp <- est <- numeric(num_MC)
  for (mc in seq_len(num_MC)) {
    dat <- make_data(base_seed + mc - 1L)
    sc  <- TRexScreeningSelector$new(
      dat$X, dat$y, seed = -1L, verbose = FALSE,
      screen_control = do.call(trex_screen_control, screen_args),
      control        = trex_ctrl
    )
    sc$select()
    sel      <- sc$selected_indices
    fdp[mc]  <- compute_fdp(sel, dat$true_support)
    tpp[mc]  <- compute_tpp(sel, dat$true_support)
    est[mc]  <- sc$estimated_FDR
  }
  cat(sprintf("  %s -- done.  TPR=%.4f  FDR=%.4f  Est.FDR=%.4f\n",
              label, mean(tpp), mean(fdp), mean(est)))
  list(avg_fdr = mean(fdp), avg_tpr = mean(tpp), avg_est_fdr = mean(est))
}

# Internal: format one right-aligned numeric row of a sweep table.
.fmt_row <- function(name, metric, vec, name_w = 31, metric_w = 10,
                     sweep_w = 10, col_w = 10) {
  s <- paste0(formatC(name, width = -name_w, flag = "-"),
              formatC(metric, width = -metric_w, flag = "-"),
              strrep(" ", sweep_w))
  s <- paste0(s, paste(formatC(sprintf("%.4f", vec), width = col_w), collapse = ""))
  s
}

#' Save + print a sweep-variable Screen-TRex MC table (mirrors
#' save_and_print_mc_results()). Writes TXT + tidy CSV and echoes to console.
#'
#' @param num_MC       Monte Carlo trials.
#' @param out_dir      Output directory.
#' @param file_stem    File stem (no extension).
#' @param x_values     Sweep grid (SNR or rho values).
#' @param method_names Ordered method names.
#' @param fdr_map,tpr_map,est_fdr_map Named lists: method -> numeric vector.
#' @param sweep_label  Column label (default "SNR").
save_and_print_scr_mc <- function(num_MC, out_dir, file_stem, x_values,
                                  method_names, fdr_map, tpr_map, est_fdr_map,
                                  sweep_label = "SNR") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- character(0)
  add <- function(...) lines[[length(lines) + 1L]] <<- paste0(...)

  add("")
  add(strrep("=", 70L))
  add(sprintf("=== Screen-TRex Results (averaged over %d Monte Carlo runs) ===", num_MC))
  add(strrep("=", 70L))
  add("")

  hdr <- paste0(formatC("Method", width = -31, flag = "-"),
                formatC("Metric", width = -10, flag = "-"),
                formatC(sweep_label, width = -10, flag = "-"),
                paste(formatC(sprintf("%.2f", x_values), width = 10), collapse = ""))
  add(hdr)
  add(strrep("-", nchar(hdr)))

  for (nm in method_names) {
    add(.fmt_row(nm, "FDR",      fdr_map[[nm]]))
    add(.fmt_row("", "TPR",      tpr_map[[nm]]))
    add(.fmt_row("", "Est. FDR", est_fdr_map[[nm]]))
    add("")
  }

  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n")
  txt_path <- file.path(out_dir, paste0(file_stem, ".txt"))
  writeLines(lines, txt_path)
  cat(sprintf("[Info] TXT results saved to: %s\n", txt_path))

  # tidy long CSV
  rows <- list()
  for (nm in method_names) {
    for (i in seq_along(x_values)) {
      rows[[length(rows) + 1L]] <- data.frame(
        method = nm, metric = "FDR",      x = x_values[i], value = fdr_map[[nm]][i])
      rows[[length(rows) + 1L]] <- data.frame(
        method = nm, metric = "TPR",      x = x_values[i], value = tpr_map[[nm]][i])
      rows[[length(rows) + 1L]] <- data.frame(
        method = nm, metric = "Est. FDR", x = x_values[i], value = est_fdr_map[[nm]][i])
    }
  }
  csv_df <- do.call(rbind, rows)
  names(csv_df)[3L] <- sweep_label
  csv_path <- file.path(out_dir, paste0(file_stem, ".csv"))
  utils::write.csv(csv_df, csv_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("[Info] CSV results saved to: %s\n\n", csv_path))
}


# ==============================================================================
# Section 5 — Biobank (Algorithm 1) helpers
# ==============================================================================

# Canonical method-name keys returned in the `method_used` column.
BIOBANK_METHODS <- c("Screen-TRex (ordinary)",
                     "Screen-TRex (bootstrap-CI)",
                     "T-Rex (fallback)")

#' Print a single-phenotype biobank result (mirrors print_biobank_single_result()).
#'
#' @param rec          A single biobank record: the return value of
#'                     selector$select() for a vector y (one BiobankScreenTRex
#'                     record, with fields method_used, estimated_FDR,
#'                     selected_indices, ...).
#' @param true_support 1-based ground-truth active indices.
.print_biobank_single <- function(rec, true_support) {
  sel <- rec$selected_indices
  cat(sprintf("\n  Method used:       %s\n", rec$method_used))
  cat(sprintf("  Estimated FDR:     %.4f\n", rec$estimated_FDR))
  cat(sprintf("  alpha_hat (ord):   %.4f\n", rec$estimated_FDR_screen_ordinary))
  cat(sprintf("  alpha_hat_C (boot):%.4f\n", rec$estimated_FDR_screen_bootstrap))
  cat(sprintf("  Selected (%d):  %s\n", length(sel), paste(sort(sel), collapse = " ")))
  cat(sprintf("  FDP: %.4f  |  TPP: %.4f\n",
              compute_fdp(sel, true_support), compute_tpp(sel, true_support)))
  cat(sprintf("  True support:   %s\n\n", paste(sort(true_support), collapse = " ")))
}

#' Print a per-phenotype biobank table (mirrors print_biobank_phenotype_table()).
#'
#' @param res           A list of biobank records (selector$select() for a
#'                      matrix Y), one record per phenotype.
#' @param true_supports Ground-truth supports (one entry per phenotype).
.print_biobank_table <- function(res, true_supports) {
  cat("\n", strrep("=", 80L), "\n", sep = "")
  cat("=== Results by Phenotype ===\n")
  cat(strrep("=", 80L), "\n\n", sep = "")
  cat(sprintf("%6s  %-26s%10s%10s%10s%10s\n",
              "Pheno", "Method", "Est.FDR", "Selected", "FDP", "TPP"))
  cat(strrep("-", 80L), "\n", sep = "")
  for (i in seq_along(true_supports)) {
    rec <- res[[i]]
    sel <- rec$selected_indices
    cat(sprintf("%6d  %-26s%10.4f%10d%10.4f%10.4f\n",
                i, rec$method_used, rec$estimated_FDR, length(sel),
                compute_fdp(sel, true_supports[[i]]),
                compute_tpp(sel, true_supports[[i]])))
  }
  cat("\n")
}

#' Print a multi-phenotype biobank summary (mirrors print_biobank_summary()).
#'
#' @param res A list of biobank records (selector$select() for a matrix Y).
.print_biobank_summary <- function(res) {
  method_used <- vapply(res, function(r) r$method_used, character(1))
  est_fdr     <- vapply(res, function(r) r$estimated_FDR, numeric(1))
  n_sel       <- vapply(res, function(r) length(r$selected_indices), integer(1))
  cat(strrep("=", 70L), "\n", sep = "")
  cat("=== Summary ===\n")
  cat(strrep("=", 70L), "\n", sep = "")
  cat(sprintf("  Total phenotypes:        %d\n", length(res)))
  cat(sprintf("  Screen-TRex (ordinary):  %d\n",
              sum(method_used == "Screen-TRex (ordinary)")))
  cat(sprintf("  Screen-TRex (bootstrap): %d\n",
              sum(method_used == "Screen-TRex (bootstrap-CI)")))
  cat(sprintf("  T-Rex (fallback):        %d\n",
              sum(method_used == "T-Rex (fallback)")))
  cat(sprintf("  Avg estimated FDR:       %.4f\n", mean(est_fdr)))
  cat(sprintf("  Avg selected count:      %.2f\n\n", mean(n_sel)))
}

#' Save + print a Biobank MC table with per-method Usage (%) rows
#' (mirrors save_and_print_biobank_mc_results()).
save_and_print_biobank_mc <- function(num_MC, out_dir, file_stem, snr_values,
                                      method_names, fdp_map, tpp_map,
                                      est_fdr_map, usage_map) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- character(0)
  add <- function(...) lines[[length(lines) + 1L]] <<- paste0(...)

  add("")
  add(strrep("=", 70L))
  add(sprintf("=== Biobank Screen-TRex Results (averaged over %d Monte Carlo runs) ===",
              num_MC))
  add(strrep("=", 70L))
  add("")

  hdr <- paste0(formatC("Method", width = -31, flag = "-"),
                formatC("Metric", width = -12, flag = "-"),
                formatC("SNR", width = -10, flag = "-"),
                paste(formatC(sprintf("%.2f", snr_values), width = 10), collapse = ""))
  add(hdr)
  add(strrep("-", nchar(hdr)))

  row <- function(name, metric, vec, dp) {
    paste0(formatC(name, width = -31, flag = "-"),
           formatC(metric, width = -12, flag = "-"),
           strrep(" ", 10L),
           paste(formatC(sprintf(paste0("%.", dp, "f"), vec), width = 10),
                 collapse = ""))
  }
  for (nm in method_names) {
    add(row(nm, "Usage (%)", usage_map[[nm]] * 100, 1))
    add(row("", "FDR",       fdp_map[[nm]],         4))
    add(row("", "TPR",       tpp_map[[nm]],         4))
    add(row("", "Est. FDR",  est_fdr_map[[nm]],     4))
    add("")
  }

  cat(paste(lines, collapse = "\n"), "\n")
  txt_path <- file.path(out_dir, paste0(file_stem, ".txt"))
  writeLines(lines, txt_path)
  cat(sprintf("[Info] TXT results saved to: %s\n", txt_path))

  rows <- list()
  for (nm in method_names) {
    for (i in seq_along(snr_values)) {
      rows[[length(rows) + 1L]] <- data.frame(method = nm, metric = "Usage",
        snr = snr_values[i], value = usage_map[[nm]][i])
      rows[[length(rows) + 1L]] <- data.frame(method = nm, metric = "FDR",
        snr = snr_values[i], value = fdp_map[[nm]][i])
      rows[[length(rows) + 1L]] <- data.frame(method = nm, metric = "TPR",
        snr = snr_values[i], value = tpp_map[[nm]][i])
      rows[[length(rows) + 1L]] <- data.frame(method = nm, metric = "Est. FDR",
        snr = snr_values[i], value = est_fdr_map[[nm]][i])
    }
  }
  csv_path <- file.path(out_dir, paste0(file_stem, ".csv"))
  utils::write.csv(do.call(rbind, rows), csv_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("[Info] CSV results saved to: %s\n\n", csv_path))
}

# ==============================================================================
# End of trex_scr_sim_utils.R
# ==============================================================================
