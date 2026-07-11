# ==============================================================================
# demo_trex_scr_05_biobank_mmap.R
# ==============================================================================
#
# Biobank Screen-TRex Selector ("Algorithm 1"), memory-mapped variant.
#
# Statistically identical to demo 04, but the design matrix X and the internal
# dummy matrices are backed by memory-mapped files (use_memory_mapping = TRUE),
# so a biobank-scale X never has to reside fully in RAM. Everything else about
# Algorithm 1 — the Ordinary -> Bootstrap-CI -> T-Rex-fallback routing and its
# per-method Usage % — is unchanged.
#
# See demo 04 for the full description of the result accessors and method-used
# labels (they are identical here).
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_05_biobank_screen_trex_mmap/
#     demo_trex_scr_05_biobank_screen_trex_mmap.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part A: single phenotype, in-memory X + use_memory_mapping = TRUE (D-mmap).
#  Part B: single + multiple phenotypes, fully memory-mapped X (temp mmap file).
#  Part C: MC SNR sweep (single phenotype) with the mmap pipeline enabled.
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

N <- 300L; P <- 1000L
NUM_MC <- 30L   # C++ uses 500; reduced for a ~1-minute runtime.

biobank_ctrl <- function() trex_biobank_control(lower_bound_FDR = 0.05,
                                                upper_bound_FDR = 0.15)

#' Run a memory-mapped biobank selector on X, Y.
run_biobank_mmap <- function(X, Y, seed = 42L, verbose = FALSE, R_boot = 1000L) {
  sel <- TRexBiobankScreeningSelector$new(
    X, Y, tFDR = 0.10, seed = seed, verbose = verbose,
    biobank_control = biobank_ctrl(),
    screen_control  = trex_screen_control(R_boot = R_boot),
    control         = make_scr_trex_ctrl(use_memory_mapping = TRUE)
  )
  sel$select()
}

#' RAII-like temp-mmap lifecycle helper (see demo 02).
with_temp_mmap <- function(X_matrix, fn) {
  x_path <- tempfile(fileext = ".dat")
  on.exit(unlink(x_path), add = TRUE)
  fn(TRexSelectorNeo::convert_to_memory_mapped(X_matrix, x_path))
}


# ==============================================================================
# Part A: single phenotype, in-memory X + D-mmap
# ==============================================================================

part_a <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part A: Biobank Screen-TRex — single phenotype, in-memory X + D-mmap\n")
  cat(strrep("=", 70L), "\n", sep = "")
  true_support <- c(5L, 28L, 43L, 150L, 399L)   # C++ 0-based {4,27,42,149,398}
  dat <- dgp_iid_snr(N, P, true_support, rep(1, length(true_support)), snr = 1.0, seed = 123L)
  res <- run_biobank_mmap(dat$X, dat$y)   # vector y -> single record
  .print_biobank_single(res, dat$true_support)
}


# ==============================================================================
# Part B: fully memory-mapped X (temp mmap file)
# ==============================================================================

part_b <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part B: Biobank Screen-TRex — memory-mapped X (single + multi phenotype)\n")
  cat(strrep("=", 70L), "\n", sep = "")

  # -- single phenotype on mmap X --
  true_support <- c(5L, 28L, 43L, 150L, 399L)
  dat <- dgp_iid_snr(N, P, true_support, rep(1, length(true_support)), snr = 1.0, seed = 123L)
  cat("Single phenotype on memory-mapped X ...\n")
  with_temp_mmap(dat$X, function(X_mmap) {
    res <- run_biobank_mmap(X_mmap, dat$y)   # vector y -> single record
    .print_biobank_single(res, dat$true_support)
  })

  # -- multiple phenotypes sharing one mmap X --
  q <- 5L
  snr_values <- c(1.0, 2.0, 5.0, 10.0, 20.0)
  sup_sizes  <- c(5L, 6L, 5L, 7L, 5L)
  set.seed(4212L)
  X <- matrix(rnorm(N * P), N, P)
  true_supports <- vector("list", q); Y <- matrix(0, N, q)
  set.seed(777L)
  for (k in seq_len(q)) {
    sup <- sample.int(P, sup_sizes[k]); true_supports[[k]] <- sort(sup)
    beta <- numeric(P); beta[sup] <- 1
    Y[, k] <- .make_response(X, beta, snr_values[k])
  }
  cat(sprintf("\n%d phenotypes on one shared memory-mapped X ...\n", q))
  with_temp_mmap(X, function(X_mmap) {
    res <- run_biobank_mmap(X_mmap, Y, R_boot = 500L)
    .print_biobank_table(res, true_supports)
    .print_biobank_summary(res)
  })
}


# ==============================================================================
# Part C: MC SNR sweep — single phenotype, mmap enabled
# ==============================================================================

part_c <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part C: Biobank Screen-TRex — MC SNR sweep (single phenotype), mmap enabled\n")
  cat(sprintf("n = %d, p = %d, s = 10, num_MC = %d\n", N, P, NUM_MC))
  cat(strrep("=", 70L), "\n", sep = "")
  support_size <- 10L; tFDR <- 0.10
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  ns <- length(snr_values)

  z <- function() setNames(lapply(BIOBANK_METHODS, function(x) numeric(ns)), BIOBANK_METHODS)
  fdp <- z(); tpp <- z(); est <- z(); usage <- z()

  for (si in seq_along(snr_values)) {
    snr <- snr_values[si]
    for (mc in seq_len(NUM_MC)) {
      seed <- 1000L * (si - 1L) + mc
      set.seed(seed + 500000L)
      sup <- sample.int(P, support_size)
      dat <- dgp_iid_snr(N, P, sup, rep(1, support_size), snr, seed = seed)
      res <- with_temp_mmap(dat$X, function(X_mmap)
        run_biobank_mmap(X_mmap, dat$y, seed = -1L, R_boot = 500L))  # single record
      sel <- res$selected_indices
      m <- res$method_used
      fdp[[m]][si]   <- fdp[[m]][si]   + compute_fdp(sel, dat$true_support)
      tpp[[m]][si]   <- tpp[[m]][si]   + compute_tpp(sel, dat$true_support)
      usage[[m]][si] <- usage[[m]][si] + 1
      est[["Screen-TRex (ordinary)"]][si]    <- est[["Screen-TRex (ordinary)"]][si]    + res$estimated_FDR_screen_ordinary
      est[["Screen-TRex (bootstrap-CI)"]][si] <- est[["Screen-TRex (bootstrap-CI)"]][si] + res$estimated_FDR_screen_bootstrap
      est[["T-Rex (fallback)"]][si]          <- est[["T-Rex (fallback)"]][si]          + tFDR
    }
    cat(sprintf("  SNR %.2f -- completed %d runs.\n", snr, NUM_MC))
  }

  for (m in BIOBANK_METHODS) {
    for (i in seq_len(ns)) {
      if (usage[[m]][i] > 0) { fdp[[m]][i] <- fdp[[m]][i] / usage[[m]][i]
                               tpp[[m]][i] <- tpp[[m]][i] / usage[[m]][i] }
      est[[m]][i]   <- est[[m]][i]   / NUM_MC
      usage[[m]][i] <- usage[[m]][i] / NUM_MC
    }
  }
  save_and_print_biobank_mc(NUM_MC, OUT_DIR,
    "d05_biobank_screen_trex_mc_snr_n300_p1000_s10_mmap",
    snr_values, BIOBANK_METHODS, fdp, tpp, est, usage)
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_a) part_a()
if (run_part_b) part_b()
if (run_part_c) part_c()
cat("\ndemo_trex_scr_05_biobank_mmap.R complete.\n")
