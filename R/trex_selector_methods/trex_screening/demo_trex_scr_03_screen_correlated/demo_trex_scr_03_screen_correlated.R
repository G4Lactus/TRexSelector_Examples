# ==============================================================================
# demo_trex_scr_03_screen_correlated.R
# ==============================================================================
#
# Screen-TRex Selector, part 03: correlated designs & dependency-aware screening.
#
# Ordinary screening thresholds each variable's voting proportion Phi_j > 0.5
# independently. Under correlated predictors, redundant correlated variables
# *split* the voting evidence, so Phi drops below 0.5 for all of them and power
# collapses. The dependency-aware (DA) variants pre-group strongly correlated
# variables before voting so the evidence is not split:
#   * TREX_DA_AR1        — for AR(1) (banded) correlation,
#   * TREX_DA_EQUI       — for equicorrelation (single shared factor),
#   * TREX_DA_BLOCK_EQUI — for block-equicorrelation (per-block factors).
# Each DA variant also reports an Estimated rho (estimated_correlation).
#
# NOTE: R's trex_screen_control() exposes only use_bootstrap_CI, R_boot,
# ci_grid_step, and trex_method. The C++ make_screen_control() additionally
# takes rho_thr_DA and n_blocks; in R those use the library's internal defaults
# (so the block DGP's n_blocks is not forwarded to the screener).
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_screen_trex_correlated/
#     demo_trex_scr_03_screen_trex_correlated.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part 1: Ordinary Screen-TRex on AR(1) data (baseline, no DA correction).
#  Part 2: DA-AR1 Screen-TRex on AR(1) data.
#  Part 3: DA-EQUI Screen-TRex on equicorrelated data.
#  Part 4: DA-BLOCK-EQUI Screen-TRex on block-equicorrelated data.
#  Part 5: MC SNR sweep on AR(1): Ordinary / Bootstrap / +DA-AR1 variants.
#  Part 6: MC SNR sweep on equicorrelated: Ordinary / Bootstrap / +DA-EQUI.
#  Part 7: MC rho sweep on AR(1) (fixed SNR): the four methods above.
#  Part 8: MC rho sweep on equicorrelated (fixed SNR).
#  Part 9: MC rho sweep on block-equicorrelated (fixed SNR).
#
# The C++ demo uses num_MC = 500; here we use a smaller count (noted below).
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================
# Setup
# ==============================================================================

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)
source(file.path(this_dir_, "..", "trex_scr_sim_utils.R"))

OUT_DIR <- file.path(this_dir_, "simulation_results", "data")

run_single_parts <- TRUE   # Parts 1-4
run_mc_parts      <- TRUE   # Parts 5-9

N <- 300L; P <- 1000L; P1 <- 10L
SNR_GRID <- c(0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
RHO_GRID <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
NUM_MC   <- 10L   # C++ uses 500; reduced so all 5 MC parts finish in ~2 min.


#' Single-run helper: build data, run Screen-TRex, print block.
run_single <- function(title, dat, trex_method, bootstrap = FALSE) {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat(title, "\n"); cat(strrep("=", 70L), "\n\n", sep = "")
  sel <- TRexScreeningSelector$new(
    dat$X, dat$y, seed = 42L, verbose = FALSE,
    screen_control = trex_screen_control(trex_method = trex_method,
                                         use_bootstrap_CI = bootstrap),
    control        = make_scr_trex_ctrl()
  )
  sel$select()
  .print_scr_result(sel, dat$true_support)
}

#' Generic MC sweep over one variable (SNR or rho), for a list of methods.
run_mc_sweep <- function(part_title, dgp_fn, sweep_name, sweep_values,
                         methods, file_stem, fixed) {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat(part_title, "\n"); cat(strrep("=", 70L), "\n\n", sep = "")
  method_names <- vapply(methods, function(m) m$name, character(1))
  init <- function() setNames(
    lapply(method_names, function(x) numeric(length(sweep_values))), method_names)
  fdr_map <- init(); tpr_map <- init(); est_fdr_map <- init()
  trex_ctrl <- make_scr_trex_ctrl()

  for (m in methods) {
    cat(strrep("=", 70L), "\n", sep = "")
    cat(sprintf("Method: %s\n", m$name)); cat(strrep("=", 70L), "\n", sep = "")
    for (vi in seq_along(sweep_values)) {
      v <- sweep_values[vi]
      make_data <- function(seed) dgp_fn(v, seed)   # v = SNR or rho
      res <- run_mc_screen(
        NUM_MC, sprintf("%s  %s=%.2f", m$name, sweep_name, v), make_data, trex_ctrl,
        screen_args = list(trex_method = m$trex_method, use_bootstrap_CI = m$bootstrap),
        base_seed = 42L + (vi - 1L) * 1000L)
      fdr_map[[m$name]][vi]     <- res$avg_fdr
      tpr_map[[m$name]][vi]     <- res$avg_tpr
      est_fdr_map[[m$name]][vi] <- res$avg_est_fdr
    }
  }
  save_and_print_scr_mc(NUM_MC, OUT_DIR, file_stem, sweep_values,
                        method_names, fdr_map, tpr_map, est_fdr_map,
                        sweep_label = sweep_name)
}


# ==============================================================================
# Parts 1-4: single runs on correlated designs
# ==============================================================================

single_parts <- function() {
  # Part 1: Ordinary on AR(1) (rho=0.7, SNR=5, beta=3) — no DA correction.
  run_single("Part 1: Ordinary Screen-TRex on AR(1) data (no DA correction)",
             dgp_ar1(N, P, P1, rho = 0.7, snr = 5.0, beta_val = 3.0, seed = 1L),
             "TREX")
  # Part 2: DA-AR1 on AR(1) (rho=0.5, SNR=1, beta=1).
  run_single("Part 2: DA-AR1 Screen-TRex on AR(1) data",
             dgp_ar1(N, P, P1, rho = 0.5, snr = 1.0, beta_val = 1.0, seed = 1L),
             "TREX_DA_AR1")
  # Part 3: DA-EQUI on equicorrelated (rho=0.4, SNR=1, beta=1).
  run_single("Part 3: DA-EQUI Screen-TRex on equicorrelated data",
             dgp_equi(N, P, P1, rho = 0.4, snr = 1.0, beta_val = 1.0, seed = 1L),
             "TREX_DA_EQUI")
  # Part 4: DA-BLOCK-EQUI on block-equicorrelated (rho=0.5, SNR=1, beta=1).
  run_single("Part 4: DA-BLOCK-EQUI Screen-TRex on block-equicorrelated data",
             dgp_block_equi(N, P, P1, rho = 0.5, n_blocks = 5L, snr = 1.0,
                            beta_val = 1.0, seed = 1L),
             "TREX_DA_BLOCK_EQUI")
}


# ==============================================================================
# Parts 5-9: Monte Carlo sweeps
# ==============================================================================

mc_parts <- function() {
  methods_ar1 <- list(
    list(name = "Screen-TRex Ordinary",         trex_method = "TREX",        bootstrap = FALSE),
    list(name = "Screen-TRex Bootstrap",        trex_method = "TREX",        bootstrap = TRUE),
    list(name = "Screen-TRex Ord+DA-AR1",       trex_method = "TREX_DA_AR1", bootstrap = FALSE),
    list(name = "Screen-TRex Boot+DA-AR1",      trex_method = "TREX_DA_AR1", bootstrap = TRUE))
  methods_equi <- list(
    list(name = "Screen-TRex Ordinary",         trex_method = "TREX",         bootstrap = FALSE),
    list(name = "Screen-TRex Bootstrap",        trex_method = "TREX",         bootstrap = TRUE),
    list(name = "Screen-TRex Ord+DA-EQUI",      trex_method = "TREX_DA_EQUI", bootstrap = FALSE),
    list(name = "Screen-TRex Boot+DA-EQUI",     trex_method = "TREX_DA_EQUI", bootstrap = TRUE))
  methods_block <- list(
    list(name = "Screen-TRex Ordinary",         trex_method = "TREX",               bootstrap = FALSE),
    list(name = "Screen-TRex Bootstrap",        trex_method = "TREX",               bootstrap = TRUE),
    list(name = "Screen-TRex Ord+DA-BLOCK",     trex_method = "TREX_DA_BLOCK_EQUI", bootstrap = FALSE),
    list(name = "Screen-TRex Boot+DA-BLOCK",    trex_method = "TREX_DA_BLOCK_EQUI", bootstrap = TRUE))

  # Part 5: SNR sweep on AR(1) (rho=0.5, beta=1)
  run_mc_sweep("Part 5: MC SNR sweep on AR(1) data (rho=0.5)",
               function(snr, seed) dgp_ar1(N, P, P1, 0.5, snr, 1.0, seed),
               "SNR", SNR_GRID, methods_ar1, "d03_da_ar1_mc_n300_p1000_rho0.50")

  # Part 6: SNR sweep on equicorrelated (rho=0.4, beta=3)
  run_mc_sweep("Part 6: MC SNR sweep on equicorrelated data (rho=0.4)",
               function(snr, seed) dgp_equi(N, P, P1, 0.4, snr, 3.0, seed),
               "SNR", SNR_GRID, methods_equi, "d03_da_equi_mc_n300_p1000_rho0.40")

  # Part 7: rho sweep on AR(1) (SNR=1, beta=1)
  run_mc_sweep("Part 7: MC rho sweep on AR(1) data (SNR=1)",
               function(rho, seed) dgp_ar1(N, P, P1, rho, 1.0, 1.0, seed),
               "rho", RHO_GRID, methods_ar1, "d03_da_ar1_rho_sweep_n300_p1000_snr1.00")

  # Part 8: rho sweep on equicorrelated (SNR=1, beta=3)
  run_mc_sweep("Part 8: MC rho sweep on equicorrelated data (SNR=1)",
               function(rho, seed) dgp_equi(N, P, P1, rho, 1.0, 3.0, seed),
               "rho", RHO_GRID, methods_equi, "d03_da_equi_rho_sweep_n300_p1000_snr1.00")

  # Part 9: rho sweep on block-equicorrelated (SNR=1, beta=3, 5 blocks)
  run_mc_sweep("Part 9: MC rho sweep on block-equicorrelated data (SNR=1)",
               function(rho, seed) dgp_block_equi(N, P, P1, rho, 5L, 1.0, 3.0, seed),
               "rho", RHO_GRID, methods_block, "d03_da_block_equi_rho_sweep_n300_p1000_snr1.00")
}


# ==============================================================================
# Main
# ==============================================================================

if (run_single_parts) single_parts()
if (run_mc_parts)      mc_parts()
cat("\ndemo_trex_scr_03_screen_correlated.R complete.\n")
