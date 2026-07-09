# ==============================================================================
# demo_trex_scr_06_screen_solvers.R
# ==============================================================================
#
# Screen-TRex Selector, part 06: solver-backend comparison.
#
# Screening thresholds the dummy-based voting proportion Phi_j, but Phi itself
# is produced by an underlying T-Rex solver. This demo asks whether the choice
# of solver backend changes screening performance:
#   * TLARS    — T-Rex Least-Angle Regression Solver (the default),
#   * TAFS     — T-Rex Adaptive Forward Selection (rho_afs = 0.3),
#   * TOMP     — T-Rex Orthogonal Matching Pursuit,
# each combined with Ordinary and Bootstrap-CI screening.
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_screen_trex_solvers/
#     demo_trex_scr_06_screen_trex_solvers.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part 1: MC SNR sweep — three solvers, Ordinary screening only (3 series).
#  Part 2: MC SNR sweep — three solvers x {Ordinary, Bootstrap-CI} (6 series).
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

OUT_DIR <- file.path(this_dir_, "simulation_results")

run_part_1 <- TRUE
run_part_2 <- TRUE

N <- 300L; P <- 1000L; SUPPORT_SIZE <- 10L
SNR_GRID <- c(0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
NUM_MC   <- 20L   # C++ uses 500; reduced for a ~1-2 minute runtime.

#' A solver x screening-variant descriptor.
sv <- function(name, solver, rho_afs = 0, bootstrap = FALSE)
  list(name = name, solver = solver, rho_afs = rho_afs, bootstrap = bootstrap)

#' MC sweep over a list of solver x method variants (SNR grid).
run_solver_sweep <- function(title, variants, file_stem) {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat(title, "\n")
  cat(sprintf("High-dimensional: n = %d, p = %d, s = %d, num_MC = %d\n",
              N, P, SUPPORT_SIZE, NUM_MC))
  cat(strrep("=", 70L), "\n\n", sep = "")

  method_names <- vapply(variants, function(v) v$name, character(1))
  init <- function() setNames(
    lapply(method_names, function(x) numeric(length(SNR_GRID))), method_names)
  fdr_map <- init(); tpr_map <- init(); est_fdr_map <- init()

  for (v in variants) {
    cat(strrep("=", 70L), "\n", sep = "")
    cat(sprintf("Variant: %s\n", v$name)); cat(strrep("=", 70L), "\n", sep = "")
    # TAFS/TOMP: enable stagnation stop as in the C++ make_trex_ctrl_for().
    trex_ctrl <- make_scr_trex_ctrl(solver = v$solver, rho_afs = v$rho_afs,
                                    tloop_stagnation_stop = TRUE,
                                    tloop_max_stagnant_steps = 5L)
    for (si in seq_along(SNR_GRID)) {
      snr <- SNR_GRID[si]
      make_data <- function(seed) {
        set.seed(seed + 500000L)
        sup <- sample.int(P, SUPPORT_SIZE)
        dgp_iid_snr(N, P, sup, rep(1, SUPPORT_SIZE), snr, seed = seed)
      }
      res <- run_mc_screen(
        NUM_MC, sprintf("%s  SNR=%.2f", v$name, snr), make_data, trex_ctrl,
        screen_args = list(trex_method = "TREX", use_bootstrap_CI = v$bootstrap),
        base_seed = 24L + (si - 1L) * 1000L)
      fdr_map[[v$name]][si]     <- res$avg_fdr
      tpr_map[[v$name]][si]     <- res$avg_tpr
      est_fdr_map[[v$name]][si] <- res$avg_est_fdr
    }
  }
  save_and_print_scr_mc(NUM_MC, OUT_DIR, file_stem, SNR_GRID,
                        method_names, fdr_map, tpr_map, est_fdr_map)
}


# ==============================================================================
# Part 1: three solvers, Ordinary screening only
# ==============================================================================

part_1 <- function() {
  variants <- list(
    sv("TLARS (Ordinary)",    "TLARS"),
    sv("TAFS-0.3 (Ordinary)", "TAFS", rho_afs = 0.3),
    sv("TOMP (Ordinary)",     "TOMP"))
  run_solver_sweep("Part 1: Solver comparison — Ordinary Screen-TRex (3 series)",
                   variants, "d06_screen_trex_solvers_mc_snr_n300_p1000_s10")
}


# ==============================================================================
# Part 2: three solvers x {Ordinary, Bootstrap-CI}
# ==============================================================================

part_2 <- function() {
  variants <- list(
    sv("TLARS (Ordinary)",     "TLARS"),
    sv("TAFS-0.3 (Ordinary)",  "TAFS", rho_afs = 0.3),
    sv("TOMP (Ordinary)",      "TOMP"),
    sv("TLARS (Bootstrap)",    "TLARS", bootstrap = TRUE),
    sv("TAFS-0.3 (Bootstrap)", "TAFS", rho_afs = 0.3, bootstrap = TRUE),
    sv("TOMP (Bootstrap)",     "TOMP", bootstrap = TRUE))
  run_solver_sweep("Part 2: Solver x method comparison (6 series)",
                   variants, "d06_screen_trex_solver_method_mc_snr_n300_p1000_s10")
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_1) part_1()
if (run_part_2) part_2()
cat("\ndemo_trex_scr_06_screen_solvers.R complete.\n")
