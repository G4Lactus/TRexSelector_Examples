# =================================================================================
# demo_trex_gvs_08.R
# =================================================================================
#
# T-Rex+GVS block-equicorrelated benchmark: HAC-discovered vs oracle groups.
#
# A 100-block design (G = 100 blocks of block_size = 5 tiling p = 500 columns),
# with 10 active blocks. Four method variants are compared over a (rho x snr)
# grid, for each of three within-block truth scenarios:
#
#   M1 = EN  + HAC-clustered groups   (groups discovered from the data)
#   M2 = EN  + oracle groups          (true block labels supplied)
#   M3 = IEN + HAC-clustered groups
#   M4 = IEN + oracle groups
#
#   Truth scenarios: INDIVIDUAL (1 active var/block), REPRESENTATIVE (2-3/block),
#                    WHOLE (all 5/block).
#
# The central comparisons are M1-vs-M2 and M3-vs-M4: how much accuracy does
# automatic group discovery (HAC) cost relative to knowing the true groups?
#
# Metrics: coordinate TPR/FDR, block_FDR, a scenario-specific block metric
# (blk_hit for INDIVIDUAL/REPRESENTATIVE, full_blk for WHOLE), T_stop, and
# M_found (number of clusters the selector formed). NOTE: the cpp benchmark also
# reports a group-"purity" diagnostic computed from the selector's internal
# per-variable cluster labels; the R6 API does not expose those labels, so
# purity is reported as 1.0 for the oracle methods (true by construction) and as
# "-" (not available) for the HAC methods — use M_found as the R-side grouping
# diagnostic instead.
#
# Console-only (mirrors cpp demo 08, which writes no result files).
#
# =================================================================================

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

source(file.path(this_dir_, "..", "trex_gvs_dgps.R"))


# ==============================================================================
# Configuration
# ==============================================================================

CFG <- list(
  n               = 200L,
  p               = 500L,
  G               = 100L,
  block_size      = 5L,
  n_active_blocks = 10L,
  tFDR            = 0.1,
  K               = 20L,
  corr_max        = 0.5,      # much lower than the 0.98 used in demos 01-07
  b               = 3.0,
  base_seed       = 2026L,
  num_MC          = 500L,
  en_solver       = "TENET_AUG",
  snr_grid        = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
  rho_grid        = c(0.5, 0.9)
)

run_part_1 <- TRUE    # single-run demo
run_part_2 <- TRUE    # full benchmark (3 scenarios x rho x snr grid)


# ==============================================================================
# Block-level metric helpers  (mirror the cpp benchmark)
# ==============================================================================

#' 1-based column indices of block `blk`.
.block_cols <- function(blk, block_size) {
  ((blk - 1L) * block_size + 1L):(blk * block_size)
}

#' Fraction of active blocks with at least one selected variable.
.block_hit_rate <- function(sel, active, block_size) {
  if (length(active) == 0L) return(0)
  hits <- sum(vapply(active, function(blk)
    any(.block_cols(blk, block_size) %in% sel), logical(1L)))
  hits / length(active)
}

#' Block-level FDP = (null blocks hit) / max(total blocks hit, 1).
.block_fdp <- function(sel, active, G, block_size) {
  total_hit <- 0L
  null_hit  <- 0L
  for (blk in seq_len(G)) {
    if (any(.block_cols(blk, block_size) %in% sel)) {
      total_hit <- total_hit + 1L
      if (!(blk %in% active)) null_hit <- null_hit + 1L
    }
  }
  null_hit / max(total_hit, 1L)
}

#' Fraction of active blocks fully contained in the selection.
.full_block_rate <- function(sel, active, block_size) {
  if (length(active) == 0L) return(0)
  fully <- sum(vapply(active, function(blk)
    all(.block_cols(blk, block_size) %in% sel), logical(1L)))
  fully / length(active)
}


# ==============================================================================
# One method on one dataset
# ==============================================================================

#' Run one GVS method variant on a dataset and collect its metrics.
#'
#' @param dat      A dgp_block_equicorr() result.
#' @param gvs_type "EN" or "IEN".
#' @param oracle   TRUE to supply true block labels (groups); FALSE for HAC.
#' @param cv_seed  Seed for the elastic-net CV fold assignment.
#' @param cfg      Benchmark configuration list.
.run_block_method <- function(dat, gvs_type, oracle, cv_seed, cfg) {
  gvs_args <- list(
    gvs_type       = gvs_type,
    corr_max       = cfg$corr_max,
    hc_linkage     = "Single",
    lambda2_method = "CV_1SE_CCD",
    cv_seed        = cv_seed)
  if (gvs_type == "EN") gvs_args$en_solver <- cfg$en_solver
  if (oracle)           gvs_args$groups    <- dat$prior_groups

  sel <- TRexSelectorNeo::TRexGVSSelector$new(
    dat$X, dat$y, tFDR = cfg$tFDR, seed = -1L, verbose = FALSE,
    gvs_control = do.call(TRexSelectorNeo::trex_gvs_control, gvs_args),
    control     = TRexSelectorNeo::trex_control(solver = "TLARS", K = cfg$K))
  sel$select()

  si <- sel$selected_indices
  ts <- dat$true_support
  list(
    coord_tpr  = TRexSelectorNeo::compute_tpp(si, ts),
    coord_fdp  = TRexSelectorNeo::compute_fdp(si, ts),
    block_hit  = .block_hit_rate(si, dat$active_blocks, dat$block_size),
    block_fdp  = .block_fdp(si, dat$active_blocks, dat$G, dat$block_size),
    full_block = .full_block_rate(si, dat$active_blocks, dat$block_size),
    purity     = if (oracle) 1.0 else NA_real_,   # HAC labels not exposed in R
    T_stop     = sel$T_stop,
    M_found    = sel$max_clusters
  )
}

#' Run all four method variants (M1-M4) on one dataset.
.run_block_all_methods <- function(dat, cv_seed, cfg) {
  list(
    M1 = .run_block_method(dat, "EN",  FALSE, cv_seed, cfg),
    M2 = .run_block_method(dat, "EN",  TRUE,  cv_seed, cfg),
    M3 = .run_block_method(dat, "IEN", FALSE, cv_seed, cfg),
    M4 = .run_block_method(dat, "IEN", TRUE,  cv_seed, cfg)
  )
}


# ==============================================================================
# One (rho, snr) cell: num_MC trials, aggregated
# ==============================================================================

.METRIC_NAMES <- c("coord_tpr", "coord_fdp", "block_hit", "block_fdp",
                   "full_block", "purity", "T_stop", "M_found")

.run_block_cell <- function(cfg, scenario, rho, snr, cell_base_seed) {
  one_trial <- function(mc) {
    trial_seed <- cell_base_seed + mc - 1L
    dat <- dgp_block_equicorr(
      n = cfg$n, p = cfg$p, G = cfg$G, block_size = cfg$block_size,
      n_active_blocks = cfg$n_active_blocks, rho = rho, snr = snr,
      scenario = scenario, seed = trial_seed, b = cfg$b)
    .run_block_all_methods(dat, cv_seed = trial_seed, cfg)
  }

  mc_indices <- seq_len(cfg$num_MC)
  if (.Platform$OS.type == "unix" && num_cores > 1L) {
    trials <- parallel::mclapply(mc_indices, one_trial, mc.cores = num_cores)
  } else {
    trials <- lapply(mc_indices, one_trial)
  }
  ok <- vapply(trials, is.list, logical(1L))
  trials <- trials[ok]

  agg_method <- function(mkey) {
    out <- lapply(.METRIC_NAMES, function(mt) {
      vals <- vapply(trials, function(t) as.numeric(t[[mkey]][[mt]]), numeric(1L))
      m <- mean(vals, na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    })
    names(out) <- .METRIC_NAMES
    out
  }
  list(M1 = agg_method("M1"), M2 = agg_method("M2"),
       M3 = agg_method("M3"), M4 = agg_method("M4"))
}


# ==============================================================================
# Grid table printer  (mirrors cpp print_block_grid_table)
# ==============================================================================

.print_block_grid <- function(scenario, rho_grid, snr_grid, results, cfg) {
  blk_metric_name <- if (scenario == "WHOLE") "full_blk" else "blk_hit"
  get_blk <- function(a) if (scenario == "WHOLE") a$full_block else a$block_hit
  methods <- c("M1(EN-C)", "M2(EN-O)", "M3(IE-C)", "M4(IE-O)")

  fmt <- function(x, d) {
    if (is.na(x)) sprintf("%9s", "-")
    else sprintf(sprintf("%%9.%df", d), x)
  }

  cat("\n", strrep("=", 90), "\n", sep = "")
  cat(sprintf("  Block Bench MC: %s  (MC=%d, tFDR=%.2f, n=%d, p=%d, G=%d, blk=%d)\n",
              scenario, cfg$num_MC, cfg$tFDR, cfg$n, cfg$p, cfg$G, cfg$block_size))
  cat(strrep("=", 90), "\n", sep = "")

  for (ir in seq_along(rho_grid)) {
    cat(sprintf("\n  rho = %.2f\n", rho_grid[ir]))
    cat(sprintf("  %-12s%7s%9s%9s%9s%9s%9s%9s%9s\n",
                "method", "SNR", "TPR", "FDR", "block_FDR",
                blk_metric_name, "purity", "T_stop", "M_found"))
    cat("  ", strrep("-", 12 + 7 + 7 * 9), "\n", sep = "")
    for (is in seq_along(snr_grid)) {
      cell <- results[[ir]][[is]]
      for (mi in seq_along(methods)) {
        a <- cell[[mi]]
        cat(sprintf("  %-12s%7.2f%s%s%s%s%s%s%s\n",
                    methods[mi], snr_grid[is],
                    fmt(a$coord_tpr, 4), fmt(a$coord_fdp, 4),
                    fmt(a$block_fdp, 4), fmt(get_blk(a), 4),
                    fmt(a$purity, 4), fmt(a$T_stop, 1), fmt(a$M_found, 1)))
      }
      cat("\n")
    }
  }
  cat(strrep("=", 90), "\n\n", sep = "")
}


# ==============================================================================
# Benchmark driver for one truth scenario
# ==============================================================================

.run_block_scenario <- function(cfg, scenario) {
  n_rho <- length(cfg$rho_grid)
  n_snr <- length(cfg$snr_grid)

  cat("\n", strrep("=", 90), "\n", sep = "")
  cat(sprintf("Block Bench MC: %s   n=%d p=%d G=%d blk=%d active=%d tFDR=%.2f K=%d MC=%d corr_max=%.2f\n",
              scenario, cfg$n, cfg$p, cfg$G, cfg$block_size, cfg$n_active_blocks,
              cfg$tFDR, cfg$K, cfg$num_MC, cfg$corr_max))
  cat(strrep("=", 90), "\n", sep = "")

  results <- vector("list", n_rho)
  for (ir in seq_len(n_rho)) {
    results[[ir]] <- vector("list", n_snr)
    for (is in seq_len(n_snr)) {
      # Seeds staggered by grid-cell index so each cell uses a distinct band.
      cell_base_seed <- cfg$base_seed +
        ((ir - 1L) * n_snr + (is - 1L)) * cfg$num_MC
      cat(sprintf("  [%s] rho=%.2f snr=%.2f  [%d/%d]  %d MC x 4 methods ...\n",
                  scenario, cfg$rho_grid[ir], cfg$snr_grid[is],
                  (ir - 1L) * n_snr + is, n_rho * n_snr, cfg$num_MC))
      results[[ir]][[is]] <- .run_block_cell(
        cfg, scenario, cfg$rho_grid[ir], cfg$snr_grid[is], cell_base_seed)
    }
  }
  .print_block_grid(scenario, cfg$rho_grid, cfg$snr_grid, results, cfg)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Block-equicorrelated GVS  |  Part 1: Single-run\n")
  cat(sprintf("n=%d, p=%d, G=%d, block_size=%d, active blocks=%d\n",
              CFG$n, CFG$p, CFG$G, CFG$block_size, CFG$n_active_blocks))
  cat("Compares HAC-discovered vs oracle groups on one INDIVIDUAL-scenario draw.\n")
  cat(strrep("=", 70), "\n\n")

  dat_p1 <- dgp_block_equicorr(
    n = CFG$n, p = CFG$p, G = CFG$G, block_size = CFG$block_size,
    n_active_blocks = CFG$n_active_blocks, rho = 0.9, snr = 2.0,
    scenario = "INDIVIDUAL", seed = CFG$base_seed, b = CFG$b)
  cat(sprintf("[Part 1] Active blocks (IDs): {%s}\n",
              paste(dat_p1$active_blocks, collapse = ", ")))
  cat(sprintf("[Part 1] True active variables: s = %d  |  sigma_y = %.4f\n\n",
              dat_p1$s, dat_p1$sigma_y))

  labels <- c(M1 = "EN  + HAC ", M2 = "EN  + oracle",
              M3 = "IEN + HAC ", M4 = "IEN + oracle")
  specs  <- list(M1 = c("EN", FALSE), M2 = c("EN", TRUE),
                 M3 = c("IEN", FALSE), M4 = c("IEN", TRUE))
  cat(sprintf("  %-14s%8s%8s%10s%9s%9s%9s\n",
              "method", "TPR", "FDR", "block_FDR", "blk_hit", "T_stop", "M_found"))
  cat("  ", strrep("-", 67), "\n", sep = "")
  for (mk in names(specs)) {
    sp <- specs[[mk]]
    r  <- .run_block_method(dat_p1, sp[1], as.logical(sp[2]),
                            cv_seed = CFG$base_seed, cfg = CFG)
    cat(sprintf("  %-14s%8.4f%8.4f%10.4f%9.4f%9.1f%9.1f\n",
                labels[[mk]], r$coord_tpr, r$coord_fdp, r$block_fdp,
                r$block_hit, r$T_stop, r$M_found))
  }
  cat("\n")
}


# ==============================================================================
# Part 2: Full benchmark — all three truth scenarios
# ==============================================================================

if (run_part_2) {
  cat("\n", strrep("=", 90), "\n", sep = "")
  cat("Block-equicorrelated benchmark\n")
  cat("methods = {EN clustered, EN oracle, IEN clustered, IEN oracle}\n")
  cat("scenarios = {INDIVIDUAL, REPRESENTATIVE, WHOLE}\n")
  cat(strrep("=", 90), "\n", sep = "")

  for (scenario in c("INDIVIDUAL", "REPRESENTATIVE", "WHOLE")) {
    .run_block_scenario(CFG, scenario)
  }
}


cat("\nBlock benchmark GVS demo complete.\n")
