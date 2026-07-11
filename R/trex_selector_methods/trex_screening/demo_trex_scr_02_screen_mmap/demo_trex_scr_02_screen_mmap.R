# ==============================================================================
# demo_trex_scr_02_screen_mmap.R
# ==============================================================================
#
# Screen-TRex Selector, part 02: memory-mapped workflow.
#
# Statistically identical to demo 01, but the dummy matrices (and, in Part B,
# the design matrix X) are backed by memory-mapped files rather than held fully
# in RAM. This is the regime Screen-TRex is built for: p so large that a plain
# in-core T-Rex run is impractical.
#   * use_memory_mapping = TRUE in trex_control() activates the internal D-mmap
#     + solver-serialization pipeline.
#   * convert_to_memory_mapped() writes an in-memory X to a memory-mapped file
#     so X itself never has to reside fully in RAM.
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_02_screen_trex_mmap/
#     demo_trex_scr_02_screen_trex_mmap.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part A: in-memory X + use_memory_mapping = TRUE (D-mmap only),
#          Ordinary and Bootstrap-CI single runs.
#  Part B: fully memory-mapped X (written to a temp mmap file) + D-mmap,
#          single Ordinary run.
#  Part C: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI) with mmap enabled,
#          verifying the mmap pipeline reproduces demo 01's screening behavior.
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

run_part_a <- TRUE
run_part_b <- TRUE
run_part_c <- TRUE

TRUE_SUPPORT <- c(28L, 150L, 399L, 5L, 43L)   # C++ 0-based {27,149,398,4,42}
NUM_MC <- 40L   # C++ uses 500; reduced for a ~1-minute runtime.


#' RAII-like temp-mmap lifecycle helper (mirrors execute_with_temp_mmap()
#' from the classical trex demos): write X to a temp mmap file, yield it, and
#' unlink on exit even if an error is thrown.
with_temp_mmap <- function(X_matrix, fn) {
  x_path <- tempfile(fileext = ".dat")
  on.exit(unlink(x_path), add = TRUE)
  fn(TRexSelectorNeo::convert_to_memory_mapped(X_matrix, x_path))
}


# ==============================================================================
# Part A: in-memory X + use_memory_mapping = TRUE
# ==============================================================================

part_a <- function() {
  n <- 300L; p <- 1000L; snr <- 5.0
  for (bootstrap in c(FALSE, TRUE)) {
    label <- if (bootstrap) "Bootstrap-CI" else "Ordinary"
    cat("\n", strrep("=", 70L), "\n", sep = "")
    cat(sprintf("Part A: in-memory X + D-mmap  |  %s screening\n", label))
    cat(sprintf("n = %d, p = %d, s = %d, SNR = %.1f\n\n",
                n, p, length(TRUE_SUPPORT), snr))

    dat <- dgp_iid_snr(n, p, TRUE_SUPPORT, rep(1, length(TRUE_SUPPORT)), snr, seed = 58L)
    sel <- TRexScreeningSelector$new(
      dat$X, dat$y, seed = 42L, verbose = FALSE,
      screen_control = trex_screen_control(trex_method = "TREX",
                                           use_bootstrap_CI = bootstrap),
      control        = make_scr_trex_ctrl(use_memory_mapping = TRUE)
    )
    sel$select()
    .print_scr_result(sel, dat$true_support)
  }
}


# ==============================================================================
# Part B: fully memory-mapped X + use_memory_mapping = TRUE
# ==============================================================================

part_b <- function() {
  n <- 300L; p <- 1000L; snr <- 5.0
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part B: memory-mapped X (temp mmap file) + D-mmap  |  Ordinary screening\n")
  cat(sprintf("n = %d, p = %d, s = %d, SNR = %.1f\n\n",
              n, p, length(TRUE_SUPPORT), snr))

  dat <- dgp_iid_snr(n, p, TRUE_SUPPORT, rep(1, length(TRUE_SUPPORT)), snr, seed = 58L)
  cat("Converting X to a memory-mapped file and running Screen-TRex ...\n")
  with_temp_mmap(dat$X, function(X_mmap) {
    sel <- TRexScreeningSelector$new(
      X_mmap, dat$y, seed = 42L, verbose = FALSE,
      screen_control = trex_screen_control(trex_method = "TREX"),
      control        = make_scr_trex_ctrl(use_memory_mapping = TRUE)
    )
    sel$select()
    .print_scr_result(sel, dat$true_support)
  })
}


# ==============================================================================
# Part C: Monte Carlo SNR sweep with mmap enabled
# ==============================================================================

part_c <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part C: Monte Carlo SNR sweep (Ordinary vs. Bootstrap-CI), mmap enabled\n")
  cat(sprintf("n = 300, p = 1000, s = 10, num_MC = %d\n\n", NUM_MC))

  n <- 300L; p <- 1000L; support_size <- 10L
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  methods <- default_scr_methods()
  method_names <- vapply(methods, function(m) m$name, character(1))

  init <- function() setNames(
    lapply(method_names, function(x) numeric(length(snr_values))), method_names)
  fdr_map <- init(); tpr_map <- init(); est_fdr_map <- init()

  trex_ctrl <- make_scr_trex_ctrl(use_memory_mapping = TRUE)

  for (m in methods) {
    cat(strrep("=", 70L), "\n", sep = "")
    cat(sprintf("Method: %s\n", m$name)); cat(strrep("=", 70L), "\n", sep = "")
    for (si in seq_along(snr_values)) {
      snr <- snr_values[si]
      make_data <- function(seed) {
        set.seed(seed + 500000L)
        sup <- sample.int(p, support_size)
        dgp_iid_snr(n, p, sup, rep(1, support_size), snr, seed = seed)
      }
      res <- run_mc_screen(
        NUM_MC, sprintf("%s  SNR=%.2f", m$name, snr), make_data, trex_ctrl,
        screen_args = list(trex_method = m$trex_method, use_bootstrap_CI = m$bootstrap),
        base_seed = 24L + (si - 1L) * 1000L)
      fdr_map[[m$name]][si]     <- res$avg_fdr
      tpr_map[[m$name]][si]     <- res$avg_tpr
      est_fdr_map[[m$name]][si] <- res$avg_est_fdr
    }
  }
  save_and_print_scr_mc(NUM_MC, OUT_DIR, "d02_screen_trex_mmap_mc_n300_p1000",
                        snr_values, method_names, fdr_map, tpr_map, est_fdr_map)
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_a) part_a()
if (run_part_b) part_b()
if (run_part_c) part_c()
cat("\ndemo_trex_scr_02_screen_mmap.R complete.\n")
