# ==============================================================================
# simulation_utils.R
# ==============================================================================
#
# Shared utility functions for the DA-TRex simulation demos.
# Sourced by each demo file — do not run standalone.
#
# Contents:
#   plot_cormat()    — Plotly heatmap of a correlation matrix.
#   .print_table()   — Print a column-aligned MC results table.
#   .print_matrix()  — Print a labelled 2-D results matrix.
#   .run_mc()        — Generic parallel MC runner (core engine).
#   .run_mc_ar1()    — Wrapper: dgp_ar1_snr   + trex+DA+AR1.
#   .run_mc_nn()     — Wrapper: dgp_nn_snr    + trex+DA+NN.
#   .run_mc_equi()   — Wrapper: dgp_bt_snr    + trex+DA+equi (rho_within==rho_between).
#   .run_mc_bt()            — Wrapper: dgp_bt_snr            + trex+DA+BT (hierarchical blocks).
#   .print_table_multi()    — Print aligned table with multiple key columns.
#   .run_mc_ar1_block()     — Wrapper: dgp_ar1_block_snr      + trex+DA+BT.
#   .run_mc_ar1_block_white() — Wrapper: dgp_ar1_block_white_snr + trex+DA+BT (p=NULL).
#   .run_mc_ht_block()      — Wrapper: dgp_ht_block_snr       + trex+DA+BT.
#   .run_mc_ht_block_white() — Wrapper: dgp_ht_block_white_snr + trex+DA+BT (p=NULL).
#   .run_mc_nn_2d()         — 2-D (row x col) grid sweep using .run_mc_nn().
#   .run_mc_hastie()           — Wrapper: dgp_hastie_snr             + trex+GVS (EN or IEN).
#   .run_mc_scattered()        — Wrapper: dgp_scattered_grouped_snr  + trex+GVS (EN or IEN).
#   .run_mc_unequal_blocks()   — Wrapper: dgp_unequal_blocks_snr     + trex+GVS (EN or IEN).
#   .run_mc_mixed_blocks()     — Wrapper: dgp_mixed_blocks_snr       + trex+GVS (EN or IEN).
#   .run_mc_random_blocks()    — Wrapper: dgp_random_blocks_snr      + trex+GVS (EN or IEN).
#   .run_mc_neg_corr()         — Wrapper: dgp_neg_corr_snr           + trex+GVS (EN or IEN).
#   .run_mc_neg_traps()        — Wrapper: dgp_neg_traps_snr          + trex+GVS (EN or IEN).
#   .run_mc_equi_blocks()      — Wrapper: dgp_equi_blocks_snr        + trex+GVS (EN or IEN).
#   .run_mc_t3_equi_blocks()   — Wrapper: dgp_t3_equi_blocks_snr     + trex+GVS (EN or IEN).
#   .run_mc_ar1_blocks()       — Wrapper: dgp_ar1_blocks_snr         + trex+GVS (EN or IEN).
#   .run_mc_arma_blocks()      — Wrapper: dgp_arma_blocks_snr        + trex+GVS (EN or IEN).
#   .run_mc_hapgen_snr()        — Wrapper: dgp_hapgen_snr             + trex+GVS (EN or IEN).
# ==============================================================================


# ==============================================================================
# Plotting
# ==============================================================================

#' Plot a correlation matrix as a Plotly heatmap.
plot_cormat <- function(cor_matrix) {
  plotly::plot_ly(
    x          = colnames(cor_matrix),
    y          = rownames(cor_matrix),
    z          = cor_matrix,
    type       = "heatmap",
    colorscale = "RdBu",
    zmin       = -1,
    zmax       = 1
  )
}


# ==============================================================================
# Output helpers
# ==============================================================================

#' Print a column-aligned MC results table.
#'
#' @param results   List of named lists, each containing at least
#'                  `param_col`, `mean_FDP`, `mean_TPP`, `sd_FDP`, `sd_TPP`.
#' @param param_col Name of the swept parameter (string), e.g. "SNR" or "rho".
#' @param param_fmt Printf format for the parameter column, e.g. "%-6.1f".
.print_table <- function(results, param_col, param_fmt = "%-10.4f") {
  # Extract the width length from param_fmt (e.g. "%-6.1f" -> 6)
  fmt_width  <- as.integer(gsub("[^0-9]", "", strsplit(param_fmt, "\\.")[[1]][1]))

  # Ensure param_col header aligns perfectly with the param_fmt width
  header_fmt <- sprintf("%%-%ds", fmt_width)

  cat(sprintf(paste0(" ", header_fmt, " %-10s %-10s %-10s %-10s\n"),
              param_col, "mean_FDP", "mean_TPP", "sd_FDP", "sd_TPP"))
  cat(strrep("-", 58), "\n")
  for (r in results) {
    cat(sprintf(paste0(" ", param_fmt, " %-10.4f %-10.4f %-10.4f %-10.4f\n"),
                r[[param_col]], r$mean_FDP, r$mean_TPP, r$sd_FDP, r$sd_TPP))
  }
  cat("\n")
}


#' Print a column-aligned MC results table with multiple key columns.
#'
#' @param results    List of named lists, each containing entries for all key
#'   columns plus mean_FDP, mean_TPP, sd_FDP, sd_TPP.
#' @param key_cols   Character vector of key column names, e.g. c("Q", "p").
#' @param key_fmts   Named character vector of printf formats per key column,
#'   e.g. c(Q = "%-4d", p = "%-6d"). Must cover all entries in key_cols.
#' @param metric_fmt Printf format for the four metric columns. Default "%-10.4f".
.print_table_multi <- function(results, key_cols, key_fmts,
                               metric_fmt = "%-10.4f") {
  metric_cols <- c("mean_FDP", "mean_TPP", "sd_FDP", "sd_TPP")

  # Extract display width for each key column from its format string
  # e.g. "%-4d" -> 4, "%-6d" -> 6
  key_widths <- vapply(key_cols, function(k) {
    as.integer(gsub("[^0-9]", "", strsplit(key_fmts[[k]], "\\.")[[1]][1]))
  }, integer(1L))

  metric_width <- as.integer(gsub("[^0-9]", "", strsplit(metric_fmt, "\\.")[[1]][1]))

  # Build header
  hdr <- " "
  for (i in seq_along(key_cols)) {
    hdr <- paste0(hdr, sprintf(sprintf("%%-%ds ", key_widths[i]), key_cols[i]))
  }
  for (col in metric_cols) {
    hdr <- paste0(hdr, sprintf(sprintf("%%-%ds ", metric_width), col))
  }
  sep_w <- sum(key_widths + 1L) + (metric_width + 1L) * length(metric_cols) + 1L

  cat(trimws(hdr, "right"), "\n")
  cat(strrep("-", sep_w), "\n")

  for (r in results) {
    row <- " "
    for (i in seq_along(key_cols)) {
      row <- paste0(row, sprintf(key_fmts[[key_cols[i]]], r[[key_cols[i]]]), " ")
    }
    for (col in metric_cols) {
      row <- paste0(row, sprintf(metric_fmt, r[[col]]), " ")
    }
    cat(trimws(row, "right"), "\n")
  }
  cat("\n")
}


#' Print a labelled 2-D results matrix.
#'
#' @param mat   Numeric matrix with named rows and columns.
#' @param label Header label printed above the matrix.
.print_matrix <- function(mat, label) {
  cat(sprintf("\n %s:\n", label))
  col_labels <- colnames(mat)
  row_labels <- rownames(mat)

  header <- sprintf(" %-10s", "")
  for (cl in col_labels) {
    header <- paste0(header, sprintf(" %-8s", cl))
  }

  cat(header, "\n")
  cat(strrep("-", nchar(header) + 2L), "\n")

  for (i in seq_len(nrow(mat))) {
    row_str <- sprintf(" %-10s", row_labels[i])
    for (j in seq_len(ncol(mat))) {
      row_str <- paste0(row_str, sprintf(" %-8.3f", mat[i, j]))
    }
    cat(row_str, "\n")
  }
  cat("\n")
}


# ==============================================================================
# Monte Carlo core
# ==============================================================================

#' Generic parallel Monte Carlo runner.
#'
#' Runs `num_MC` independent trials in parallel. Each trial generates data via
#' `dgp_fun`, runs `TRexSelector::trex`, and records FDP / TPP.
#'
#' @param n,p,support,amplitude,snr  DGP parameters. Pass `support=NULL` to
#'   redraw a random support each trial. All DGP-specific parameters (e.g.
#'   `rho`, `kappa`, `n_blocks`) must be supplied via `dgp_extra_args`.
#' @param s               Number of active predictors; used only when
#'   `support=NULL`. If `NULL`, falls back to `MC_BASE$s` (for callers that
#'   define `MC_BASE` at top level).
#' @param tFDR,K          TRex parameters.
#' @param num_MC,seed     Number of MC trials and base seed.
#' @param dgp_fun         DGP function, e.g. `dgp_ar1_snr` or `dgp_bt_snr`.
#' @param dgp_extra_args  Named list forwarded to `dgp_fun` beyond the
#'   standard n/p/support/amplitude/snr/seed set.
#' @param trex_method     `method` argument for `TRexSelector::trex`.
#' @param trex_extra_args Named list forwarded to `TRexSelector::trex` beyond
#'   X/y/tFDR/K/method/verbose/seed.
#' @param n_cores         Number of parallel worker cores.
#'
#' @return list(mean_FDP, mean_TPP, sd_FDP, sd_TPP)
.run_mc <- function(
  n,
  p               = NULL,
  support,
  amplitude,
  s               = NULL,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  dgp_fun,
  dgp_extra_args  = list(),
  trex_method,
  trex_extra_args = list(),
  n_cores         = num_cores
) {

  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L

    trial_support <- if (is.null(support)) {
      s_val <- if (!is.null(s)) s else MC_BASE$s
      make_support_random(s_val, p, trial_seed)
    } else {
      support
    }

    base_dgp_args <- list(n = n, support = trial_support,
                          amplitude = amplitude, snr = snr, seed = trial_seed)
    if (!is.null(p)) base_dgp_args$p <- p
    dgp_args <- c(base_dgp_args, dgp_extra_args)
    dat <- do.call(dgp_fun, dgp_args)

    trex_args <- c(
      list(X = dat$X, y = dat$y, tFDR = tFDR, K = K,
           method = trex_method, verbose = FALSE, seed = trial_seed),
      trex_extra_args
    )
    res <- do.call(TRexSelector::trex, trex_args)

    c(
      FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta)
    )
  }

  cat(sprintf(" Running %d Monte Carlo trials on %d cores ...\n", num_MC, n_cores))

  if (.Platform$OS.type == "unix") {
    res_list <- parallel::mclapply(
      X        = seq_len(num_MC),
      FUN      = one_trial,
      mc.cores = n_cores
    )
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterSetRNGStream(cl, iseed = seed)
    parallel::clusterEvalQ(cl, {
      library(TRexSelector)
      NULL
    })
    parallel::clusterExport(
      cl,
      varlist = c("MC_BASE", "s", "support", "p", "n", "amplitude", "snr",
                  "tFDR", "K", "seed", "dgp_fun", "dgp_extra_args",
                  "trex_method", "trex_extra_args", "make_support_random"),
      envir = environment()
    )

    res_list <- parallel::parLapply(
      cl,
      X   = seq_len(num_MC),
      fun = one_trial
    )
  }

  res_mat <- do.call(rbind, res_list)
  fdps    <- res_mat[, "FDP"]
  tpps    <- res_mat[, "TPP"]

  cat(sprintf(" done. mean TPP=%.3f mean FDP=%.3f\n\n", mean(tpps), mean(fdps)))

  list(
    mean_FDP = mean(fdps),
    mean_TPP = mean(tpps),
    sd_FDP   = sd(fdps),
    sd_TPP   = sd(tpps)
  )
}


# ==============================================================================
# Named wrappers
# ==============================================================================

#' MC runner for the AR(1) DGP with trex+DA+AR1.
#'
#' Thin wrapper around `.run_mc()` fixing dgp_fun = dgp_ar1_snr and
#' trex_method = "trex+DA+AR1" with the standard DA options.
.run_mc_ar1 <- function(
  n,
  p,
  support,
  amplitude,
  rho,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n               = n,
    p               = p,
    support         = support,
    amplitude       = amplitude,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_ar1_snr,
    dgp_extra_args  = list(rho = rho),
    trex_method     = "trex+DA+AR1",
    trex_extra_args = list(cor_coef = NA, type = "lar", rho_thr_DA = 0.02),
    n_cores         = n_cores
  )
}


#' MC runner for the MA(kappa) / NN DGP with trex+DA+NN.
#'
#' Thin wrapper around `.run_mc()` fixing dgp_fun = dgp_nn_snr and
#' trex_method = "trex+DA+NN".
.run_mc_nn <- function(
  n,
  p,
  support,
  amplitude,
  kappa,
  rho,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n              = n,
    p              = p,
    support        = support,
    amplitude      = amplitude,
    snr            = snr,
    tFDR           = tFDR,
    K              = K,
    num_MC         = num_MC,
    seed           = seed,
    dgp_fun        = dgp_nn_snr,
    dgp_extra_args = list(rho = rho, kappa = kappa),
    trex_method    = "trex+DA+NN",
    n_cores        = n_cores
  )
}


#' MC runner for the pure equicorrelation DGP with trex+DA+equi.
#'
#' Uses `dgp_bt_snr` with `rho_within = rho_between = rho`, which reduces the
#' factor model to a single global factor — i.e. pure equicorrelation.
.run_mc_equi <- function(
  n,
  p,
  support,
  amplitude,
  s         = NULL,
  rho,
  n_blocks  = 10L,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores   = num_cores
) {
  .run_mc(
    n              = n,
    p              = p,
    support        = support,
    amplitude      = amplitude,
    s              = s,
    snr            = snr,
    tFDR           = tFDR,
    K              = K,
    num_MC         = num_MC,
    seed           = seed,
    dgp_fun        = dgp_bt_snr,
    dgp_extra_args = list(n_blocks    = n_blocks,
                          rho_within  = rho,
                          rho_between = rho),
    trex_method    = "trex+DA+equi",
    n_cores        = n_cores
  )
}


#' MC runner for the hierarchical-block (BT) DGP with trex+DA+BT.
#'
#' Uses `dgp_bt_snr` with independent `rho_within` and `rho_between`.
#' The `hc_dist` argument controls the HAC linkage passed to trex.
.run_mc_bt <- function(
  n,
  p,
  support,
  amplitude,
  s           = NULL,
  n_blocks    = 10L,
  rho_within,
  rho_between,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  hc_dist     = "single",
  n_cores     = num_cores
) {
  .run_mc(
    n               = n,
    p               = p,
    support         = support,
    amplitude       = amplitude,
    s               = s,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_bt_snr,
    dgp_extra_args  = list(n_blocks    = n_blocks,
                           rho_within  = rho_within,
                           rho_between = rho_between),
    trex_method     = "trex+DA+BT",
    trex_extra_args = list(hc_dist = hc_dist),
    n_cores         = n_cores
  )
}


#' MC runner for the pure AR(1)-block DGP with trex+DA+BT.
#'
#' Thin wrapper around `.run_mc()` fixing dgp_fun = dgp_ar1_block_snr and
#' trex_method = "trex+DA+BT".
.run_mc_ar1_block <- function(
  n,
  p,
  support,
  amplitude,
  n_blocks,
  rho,
  hc_dist,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n               = n,
    p               = p,
    support         = support,
    amplitude       = amplitude,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_ar1_block_snr,
    dgp_extra_args  = list(n_blocks = n_blocks, rho = rho),
    trex_method     = "trex+DA+BT",
    trex_extra_args = list(hc_dist = hc_dist),
    n_cores         = n_cores
  )
}


#' MC runner for the AR(1)-block + white-noise DGP with trex+DA+BT.
#'
#' p_ar and p_white are passed via dgp_extra_args; p is left NULL in .run_mc()
#' because the DGP signature uses p_ar/p_white rather than a single p.
#' support must always be provided (non-NULL) when using this wrapper.
.run_mc_ar1_block_white <- function(
  n,
  p_ar,
  p_white,
  support,
  amplitude,
  n_blocks,
  rho,
  hc_dist,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n               = n,
    p               = NULL,
    support         = support,
    amplitude       = amplitude,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_ar1_block_white_snr,
    dgp_extra_args  = list(p_ar = p_ar, p_white = p_white,
                           n_blocks = n_blocks, rho = rho),
    trex_method     = "trex+DA+BT",
    trex_extra_args = list(hc_dist = hc_dist),
    n_cores         = n_cores
  )
}


#' MC runner for the heavy-tailed AR(1)-block DGP with trex+DA+BT.
#'
#' Thin wrapper around `.run_mc()` fixing dgp_fun = dgp_ht_block_snr and
#' trex_method = "trex+DA+BT".
.run_mc_ht_block <- function(
  n,
  p,
  support,
  amplitude,
  n_blocks,
  rho,
  nu,
  heavy_noise,
  hc_dist,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n               = n,
    p               = p,
    support         = support,
    amplitude       = amplitude,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_ht_block_snr,
    dgp_extra_args  = list(n_blocks = n_blocks, rho = rho,
                           nu = nu, heavy_noise = heavy_noise),
    trex_method     = "trex+DA+BT",
    trex_extra_args = list(hc_dist = hc_dist),
    n_cores         = n_cores
  )
}


#' MC runner for the heavy-tailed AR(1)-block + white-noise DGP with trex+DA+BT.
#'
#' p_ar and p_white are passed via dgp_extra_args; p is left NULL in .run_mc()
#' because the DGP signature uses p_ar/p_white rather than a single p.
#' support must always be provided (non-NULL) when using this wrapper.
.run_mc_ht_block_white <- function(
  n,
  p_ar,
  p_white,
  support,
  amplitude,
  n_blocks,
  rho,
  nu,
  heavy_noise,
  hc_dist,
  snr,
  tFDR,
  K,
  num_MC,
  seed,
  n_cores = num_cores
) {
  .run_mc(
    n               = n,
    p               = NULL,
    support         = support,
    amplitude       = amplitude,
    snr             = snr,
    tFDR            = tFDR,
    K               = K,
    num_MC          = num_MC,
    seed            = seed,
    dgp_fun         = dgp_ht_block_white_snr,
    dgp_extra_args  = list(p_ar = p_ar, p_white = p_white,
                           n_blocks = n_blocks, rho = rho,
                           nu = nu, heavy_noise = heavy_noise),
    trex_method     = "trex+DA+BT",
    trex_extra_args = list(hc_dist = hc_dist),
    n_cores         = n_cores
  )
}


# ==============================================================================
# 2-D sweep helper (NN)
# ==============================================================================

#' Run a 2-D (row x col) grid sweep using .run_mc_nn().
#'
#' @param row_grid,col_grid Numeric/integer vectors defining the sweep axes.
#' @param row_name,col_name  Axis label strings (e.g. "kappa", "rho").
#' @param support            Fixed support vector, or NULL for per-trial random.
#' @param amplitude          Signal amplitude.
#' @param kappa_fun          Function(row_val, col_val) -> kappa.
#' @param rho_fun            Function(row_val, col_val) -> rho.
#' @param snr_fun            Function(row_val, col_val) -> snr.
#' @param tFDR,K,num_MC,seed TRex / MC parameters.
#'
#' @return list(mean_TPP = matrix, mean_FDP = matrix)
.run_mc_nn_2d <- function(
  row_grid,
  col_grid,
  row_name,
  col_name,
  support,
  amplitude,
  kappa_fun,
  rho_fun,
  snr_fun,
  tFDR,
  K,
  num_MC,
  seed
) {
  n_row <- length(row_grid)
  n_col <- length(col_grid)

  row_nms <- paste0(row_name, "=",
                    if (is.integer(row_grid)) row_grid
                    else sprintf("%.2f", row_grid))
  col_nms <- paste0(col_name, "=",
                    if (is.integer(col_grid)) col_grid
                    else sprintf("%.2f", col_grid))

  mat_TPP <- matrix(NA_real_, nrow = n_row, ncol = n_col,
                    dimnames = list(row_nms, col_nms))
  mat_FDP <- matrix(NA_real_, nrow = n_row, ncol = n_col,
                    dimnames = list(row_nms, col_nms))

  total_cells <- n_row * n_col
  cell_idx    <- 0L

  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      cell_idx <- cell_idx + 1L

      row_val   <- row_grid[i]
      col_val   <- col_grid[j]
      kappa_val <- kappa_fun(row_val, col_val)
      rho_val   <- rho_fun(row_val, col_val)
      snr_val   <- snr_fun(row_val, col_val)

      prefix <- sprintf(" [%2d/%2d] %s=%s %s=%s",
                        cell_idx, total_cells,
                        row_name, as.character(row_val),
                        col_name, as.character(col_val))
      cat(sprintf("%s running %d MC trials ...\n", prefix, num_MC))

      r <- .run_mc_nn(
        n         = MC_BASE$n,
        p         = MC_BASE$p,
        support   = support,
        amplitude = amplitude,
        kappa     = kappa_val,
        rho       = rho_val,
        snr       = snr_val,
        tFDR      = tFDR,
        K         = K,
        num_MC    = num_MC,
        seed      = seed
      )

      mat_TPP[i, j] <- r$mean_TPP
      mat_FDP[i, j] <- r$mean_FDP

      cat(sprintf(" %s done. TPP=%.3f FDP=%.3f\n\n",
                  prefix, mat_TPP[i, j], mat_FDP[i, j]))
    }
  }

  list(mean_TPP = mat_TPP, mean_FDP = mat_FDP)
}


# ==============================================================================
# GVS wrappers
# ==============================================================================

#' MC runner for the Hastie (Zou & Hastie 2005) DGP with trex+GVS.
#'
#' Each trial generates data via `dgp_hastie_snr(n, p, snr, sd_x, seed)` and
#' applies `TRexSelector::trex()` with `method = "trex+GVS"`.
#'
#' @param n,p         Data dimensions.
#' @param snr         Target signal-to-noise ratio.
#' @param sd_x        Within-group noise sd (rho = 1/(1+sd_x^2)).
#' @param tFDR        Target FDR level.
#' @param K           Number of random experiments per selector run.
#' @param num_MC      Number of Monte Carlo trials.
#' @param seed        Base RNG seed; trial k uses seed + k - 1.
#' @param GVS_type    GVS variant: "EN" or "IEN".
#' @param corr_max    HAC clustering threshold passed to trex.
#' @param hc_dist     HAC linkage type (default "single").
#' @param n_cores     Number of parallel cores (default num_cores).
#'
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_hastie <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_hastie_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                          seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_hastie_snr", "n", "p", "snr",
                                        "sd_x", "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


#' MC runner for the scattered-grouped DGP with trex+GVS.
#'
#' Each trial generates data via `dgp_scattered_grouped_snr()` and applies
#' `TRexSelector::trex()` with `method = "trex+GVS"`.
#' Group indices are re-randomised each trial (controlled by trial_seed).
#'
#' @param n,p         Data dimensions.
#' @param snr         Target signal-to-noise ratio.
#' @param sd_x        Within-group noise sd (rho = 1/(1+sd_x^2)).
#' @param n_groups    Number of active groups. Default 3.
#' @param group_size  Variables per group. Default 50.
#' @param tFDR        Target FDR level.
#' @param K           Number of random experiments per selector run.
#' @param num_MC      Number of Monte Carlo trials.
#' @param seed        Base RNG seed; trial k uses seed + k - 1.
#' @param GVS_type    GVS variant: "EN" or "IEN".
#' @param corr_max    HAC clustering threshold passed to trex.
#' @param hc_dist     HAC linkage type (default "single").
#' @param n_cores     Number of parallel cores (default num_cores).
#'
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_scattered <- function(
  n,
  p,
  snr,
  sd_x,
  n_groups   = 3L,
  group_size = 50L,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_scattered_grouped_snr(
      n          = n,
      p          = p,
      snr        = snr,
      sd_x       = sd_x,
      n_groups   = n_groups,
      group_size = group_size,
      seed       = trial_seed
    )
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_scattered_grouped_snr",
                                        "n", "p", "snr", "sd_x",
                                        "n_groups", "group_size",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_unequal_blocks
# ==============================================================================

#' MC runner for the unequal-blocks DGP with trex+GVS.
#'
#' Three contiguous blocks (sizes 20, 50, 80 = 150 active variables).
#' Each trial calls `dgp_unequal_blocks_snr()` then `TRexSelector::trex()`.
#'
#' @param n,p        Data dimensions.
#' @param snr        Target SNR.
#' @param sd_x       Within-block noise sd (rho = 1/(1+sd_x^2)).
#' @param tFDR       Target FDR level.
#' @param K          Random experiments per run.
#' @param num_MC     Monte Carlo trials.
#' @param seed       Base RNG seed.
#' @param GVS_type   "EN" or "IEN".
#' @param corr_max   HAC clustering threshold.
#' @param hc_dist    HAC linkage (default "single").
#' @param n_cores    Parallel cores (default num_cores).
#'
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_unequal_blocks <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_unequal_blocks_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                                  seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_unequal_blocks_snr",
                                        "n", "p", "snr", "sd_x",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_mixed_blocks
# ==============================================================================

#' MC runner for the mixed-blocks DGP with trex+GVS.
#'
#' 3 active blocks (20, 50, 80) + 1 inactive block (65), randomly ordered.
#' Active set s = 150.
#'
#' @inheritParams .run_mc_unequal_blocks
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_mixed_blocks <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_mixed_blocks_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                                seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_mixed_blocks_snr",
                                        "n", "p", "snr", "sd_x",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_random_blocks
# ==============================================================================

#' MC runner for the random-blocks DGP with trex+GVS.
#'
#' 3 active blocks (20, 50, 80) placed in a random order with noise gaps.
#' Active set s = 150.
#'
#' @inheritParams .run_mc_unequal_blocks
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_random_blocks <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_random_blocks_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                                 seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_random_blocks_snr",
                                        "n", "p", "snr", "sd_x",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_neg_corr
# ==============================================================================

#' MC runner for the negative-correlation DGP with trex+GVS.
#'
#' Active group 1-100 (+Z1/-Z1, beta = +3/-3); inactive trap 101-200.
#' Active set s = 100.
#'
#' @inheritParams .run_mc_unequal_blocks
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_neg_corr <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_neg_corr_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                            seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_neg_corr_snr",
                                        "n", "p", "snr", "sd_x",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_neg_traps
# ==============================================================================

#' MC runner for the negative-traps DGP with trex+GVS.
#'
#' Active group 1-100 (+Z1/-Z1, beta = +3/-3); two inactive traps (100 + 60).
#' Active set s = 100.
#'
#' @inheritParams .run_mc_unequal_blocks
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_neg_traps <- function(
  n,
  p,
  snr,
  sd_x,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_neg_traps_snr(n = n, p = p, snr = snr, sd_x = sd_x,
                             seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_neg_traps_snr",
                                        "n", "p", "snr", "sd_x",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_equi_blocks
# ==============================================================================

#' MC runner for the equi-correlated blocks DGP with trex+GVS.
#'
#' 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
#' Within-block equi-correlation = rho. Active set s = 150.
#'
#' @param n,p       Dataset dimensions.
#' @param snr       Signal-to-noise ratio.
#' @param rho       Within-block equi-correlation.
#' @param tFDR,K    TRex parameters.
#' @param num_MC,seed Monte Carlo parameters.
#' @param GVS_type  "EN" or "IEN".
#' @param corr_max  Hierarchical clustering threshold.
#' @param hc_dist   HC linkage. Default "single".
#' @param n_cores   Parallel worker count.
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_equi_blocks <- function(
  n,
  p,
  snr,
  rho,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_equi_blocks_snr(n = n, p = p, snr = snr, rho = rho,
                               seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_equi_blocks_snr",
                                        "n", "p", "snr", "rho",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_t3_equi_blocks
# ==============================================================================

#' MC runner for the heavy-tailed equi-correlated blocks DGP with trex+GVS.
#'
#' Same block structure as .run_mc_equi_blocks but with Student-t(df) noise.
#' Active set s = 150.
#'
#' @inheritParams .run_mc_equi_blocks
#' @param df Degrees of freedom for the t distribution. Default 3.
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_t3_equi_blocks <- function(
  n,
  p,
  snr,
  rho,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  df       = 3L,
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_t3_equi_blocks_snr(n = n, p = p, snr = snr, rho = rho,
                                   df = df, seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_t3_equi_blocks_snr",
                                        "n", "p", "snr", "rho", "df",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_ar1_blocks
# ==============================================================================

#' MC runner for the block-structured AR(1) DGP with trex+GVS.
#'
#' 4 blocks (sizes 20, 50, 80, 65): 3 active (beta=3), 1 inactive trap.
#' Within-block AR(1) correlation: Cor(j,k) = rho^|j-k|. Active set s = 150.
#'
#' @inheritParams .run_mc_equi_blocks
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_ar1_blocks <- function(
  n,
  p,
  snr,
  rho,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_ar1_blocks_snr(n = n, p = p, snr = snr, rho = rho,
                              seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_ar1_blocks_snr",
                                        "n", "p", "snr", "rho",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_arma_blocks
# ==============================================================================

#' MC runner for the ARMA-structured blocks DGP with trex+GVS.
#'
#' 4 blocks (AR1 / MA3 / ARMA21 / AR1-trap): 3 active (beta=3), 1 inactive.
#' Columns column-standardized. Active set s = 150.
#'
#' @param n,p       Dataset dimensions.
#' @param snr       Signal-to-noise ratio.
#' @param ar_coef   AR coefficient shared by AR(1) and ARMA(2,1) blocks.
#' @param tFDR,K    TRex parameters.
#' @param num_MC,seed Monte Carlo parameters.
#' @param GVS_type  "EN" or "IEN".
#' @param corr_max  Hierarchical clustering threshold.
#' @param hc_dist   HC linkage. Default "single".
#' @param n_cores   Parallel worker count.
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_arma_blocks <- function(
  n,
  p,
  snr,
  ar_coef,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_arma_blocks_snr(n = n, p = p, snr = snr, ar_coef = ar_coef,
                               seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_arma_blocks_snr",
                                        "n", "p", "snr", "ar_coef",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, library(TRexSelector))
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}


# ==============================================================================
# .run_mc_hapgen_snr
# ==============================================================================
#' MC runner for the quasi-hapgen LD-block DGP with trex+GVS.
#'
#' 7 irregular LD blocks (p = 500, fixed layout). Active set: blocks 1, 3, 4 (s = 130).
#' Two weak long-range cross-block LD pairs. Sigma fixed per rho_scale, reused across trials.
#'
#' @param n         Number of observations.
#' @param snr       Signal-to-noise ratio.
#' @param rho_scale Overall correlation scale in [0, 1].
#' @param tFDR,K    TRex parameters.
#' @param num_MC,seed Monte Carlo parameters.
#' @param GVS_type  "EN" or "IEN".
#' @param corr_max  Hierarchical clustering threshold.
#' @param hc_dist   HC linkage. Default "single".
#' @param n_cores   Parallel worker count.
#' @return Named list: mean_FDP, mean_TPP, sd_FDP, sd_TPP.
.run_mc_hapgen_snr <- function(
  n,
  snr,
  rho_scale,
  tFDR,
  K,
  num_MC,
  seed,
  GVS_type,
  corr_max,
  hc_dist  = "single",
  n_cores  = num_cores
) {
  one_trial <- function(mc) {
    trial_seed <- seed + mc - 1L
    dat <- dgp_hapgen_snr(n = n, snr = snr, rho_scale = rho_scale,
                          seed = trial_seed)
    res <- TRexSelector::trex(
      X             = dat$X,
      y             = dat$y,
      tFDR          = tFDR,
      K             = K,
      method        = "trex+GVS",
      GVS_type      = GVS_type,
      lambda_2_lars = NULL,
      hc_dist       = hc_dist,
      corr_max      = corr_max,
      verbose       = FALSE,
      seed          = trial_seed
    )
    c(FDP = TRexSelector::FDP(res$selected_var, dat$beta),
      TPP = TRexSelector::TPP(res$selected_var, dat$beta))
  }

  mc_indices <- seq_len(num_MC)

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    results <- parallel::mclapply(mc_indices, one_trial, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,
                            varlist = c("dgp_hapgen_snr",
                                        "n", "snr", "rho_scale",
                                        "tFDR", "K", "GVS_type",
                                        "hc_dist", "corr_max", "seed"),
                            envir   = environment())
    parallel::clusterEvalQ(cl, {
      library(TRexSelector)
      library(Matrix)
    })
    results <- parallel::parLapply(cl, mc_indices, one_trial)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(mc_indices, one_trial)
  }

  fdp_vec <- sapply(results, `[[`, "FDP")
  tpp_vec <- sapply(results, `[[`, "TPP")

  list(
    mean_FDP = mean(fdp_vec),
    mean_TPP = mean(tpp_vec),
    sd_FDP   = sd(fdp_vec),
    sd_TPP   = sd(tpp_vec)
  )
}
