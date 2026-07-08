# ==============================================================================
# demo_trex_scr_04_biobank_inmem.R
# ==============================================================================
#
# Biobank Screen-TRex Selector ("Algorithm 1"), in-memory variant.
#
# The biobank workflow screens MANY phenotypes against ONE shared design matrix
# X. For each phenotype it routes adaptively:
#   1. try Screen-TRex Ordinary; if its estimated FDR lands in the acceptable
#      window [lower_bound_FDR, upper_bound_FDR] it is accepted;
#   2. otherwise try Screen-TRex Bootstrap-CI;
#   3. otherwise fall back to the classical T-Rex selector at tFDR.
# The fraction of phenotypes routed to each method is the "Usage %".
#
# Result accessors (from selector$select(), discovered at runtime):
#   $statistics                 data.frame, one row per phenotype:
#       phenotype_index, estimated_FDR, method_used,
#       estimated_FDR_screen_ordinary, estimated_FDR_screen_bootstrap,
#       used_fallback_trex
#   $selected_indices           list, 1-based selected indices per phenotype
#   $selected_indices_screen_ordinary / $..._screen_bootstrap
#   method_used values: "Screen-TRex (ordinary)", "Screen-TRex (bootstrap-CI)",
#                       "T-Rex (fallback)".
#
# The file mirrors:
# cpp/trex_selector_methods/trex_screening/demo_trex_scr_04_biobank_screen_trex_inmem/
#     demo_trex_scr_04_biobank_screen_trex_inmem.cpp
#
# Demo content:
# ----------------------------------------------------------------
#  Part 1: single phenotype — basic Algorithm 1 functionality.
#  Part 2: multiple phenotypes (q=5) sharing one X, per-phenotype routing table.
#  Part 3: MC SNR sweep (single phenotype): per-method Usage %, FDR, TPR, Est.FDR.
#  Part 4: MC SNR sweep (multiple phenotypes): same per-method metrics.
#
# The C++ demo uses num_MC = 500; here we use smaller counts (noted below).
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
source(file.path(this_dir_, "trex_scr_sim_utils.R"))

OUT_DIR <- file.path(this_dir_, "simulation_results")

run_part_1 <- TRUE
run_part_2 <- TRUE
run_part_3 <- TRUE
run_part_4 <- TRUE

N <- 300L; P <- 1000L
NUM_MC_SINGLE <- 40L   # C++ uses 500; reduced for a ~1-minute runtime.
NUM_MC_MULTI  <- 15L   # C++ uses 500; reduced (q=5 phenotypes per run).

biobank_ctrl <- function() trex_biobank_control(lower_bound_FDR = 0.05,
                                                upper_bound_FDR = 0.15)

#' Build a biobank selector, run it, return the select() result list.
run_biobank <- function(X, Y, seed = 42L, verbose = FALSE, R_boot = 1000L) {
  sel <- TRexBiobankScreeningSelector$new(
    X, Y, tFDR = 0.10, seed = seed, verbose = verbose,
    biobank_control = biobank_ctrl(),
    screen_control  = trex_screen_control(R_boot = R_boot),
    control         = make_scr_trex_ctrl()
  )
  sel$select()
}


# ==============================================================================
# Part 1: single phenotype
# ==============================================================================

part_1 <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 1: Biobank Screen-TRex — Single Phenotype (in-memory)\n")
  cat(strrep("=", 70L), "\n", sep = "")
  # C++ 0-based support {4,27,42,149,398} -> R 1-based
  true_support <- c(5L, 28L, 43L, 150L, 399L)
  snr <- 1.0
  cat(sprintf("DGP: n = %d, p = %d, p1 = %d, SNR = %.1f\n",
              N, P, length(true_support), snr))
  dat <- dgp_iid_snr(N, P, true_support, rep(1, length(true_support)), snr, seed = 123L)
  res <- run_biobank(dat$X, matrix(dat$y, N, 1L))
  .print_biobank_single(res, dat$true_support)
}


# ==============================================================================
# Part 2: multiple phenotypes sharing one design matrix
# ==============================================================================

part_2 <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 2: Biobank Screen-TRex — Multiple Phenotypes (in-memory)\n")
  cat(strrep("=", 70L), "\n", sep = "")
  q <- 5L
  snr_values <- c(1.0, 2.0, 5.0, 10.0, 20.0)
  sup_sizes  <- c(5L, 6L, 5L, 7L, 5L)

  # Shared X (fixed seed), per-phenotype response with its own support/SNR.
  set.seed(4212L)
  X <- matrix(rnorm(N * P), N, P)
  true_supports <- vector("list", q)
  Y <- matrix(0, N, q)
  set.seed(777L)
  for (k in seq_len(q)) {
    sup <- sample.int(P, sup_sizes[k])
    true_supports[[k]] <- sort(sup)
    beta <- numeric(P); beta[sup] <- 1
    Y[, k] <- .make_response(X, beta, snr_values[k])
    cat(sprintf("Phenotype %d: p1 = %d, SNR = %.2f\n", k, sup_sizes[k], snr_values[k]))
  }
  cat(sprintf("\nRunning Algorithm 1 for %d phenotypes ...\n", q))
  res <- run_biobank(X, Y, R_boot = 500L)
  .print_biobank_table(res, true_supports)
  .print_biobank_summary(res)
}


# ==============================================================================
# MC accumulation helper
# ==============================================================================

# Accumulate one phenotype's contribution into the running maps (by reference
# via <<- in the caller's environment is avoided; we return updated maps).
.accumulate <- function(maps, st_row, sel, true_support, si, tFDR) {
  m <- st_row$method_used
  maps$fdp[[m]][si]   <- maps$fdp[[m]][si]   + compute_fdp(sel, true_support)
  maps$tpp[[m]][si]   <- maps$tpp[[m]][si]   + compute_tpp(sel, true_support)
  maps$usage[[m]][si] <- maps$usage[[m]][si] + 1
  maps$est[["Screen-TRex (ordinary)"]][si]     <-
    maps$est[["Screen-TRex (ordinary)"]][si]     + st_row$estimated_FDR_screen_ordinary
  maps$est[["Screen-TRex (bootstrap-CI)"]][si]  <-
    maps$est[["Screen-TRex (bootstrap-CI)"]][si]  + st_row$estimated_FDR_screen_bootstrap
  maps$est[["T-Rex (fallback)"]][si]           <-
    maps$est[["T-Rex (fallback)"]][si]           + tFDR
  maps
}

.new_maps <- function(ns) {
  z <- function() setNames(lapply(BIOBANK_METHODS, function(x) numeric(ns)), BIOBANK_METHODS)
  list(fdp = z(), tpp = z(), est = z(), usage = z())
}


# ==============================================================================
# Part 3: MC SNR sweep — single phenotype
# ==============================================================================

part_3 <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 3: Biobank Screen-TRex — MC SNR sweep, single phenotype (in-memory)\n")
  cat(sprintf("n = %d, p = %d, s = 10, num_MC = %d\n", N, P, NUM_MC_SINGLE))
  cat(strrep("=", 70L), "\n", sep = "")
  support_size <- 10L; tFDR <- 0.10
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  ns <- length(snr_values)
  maps <- .new_maps(ns)

  for (si in seq_along(snr_values)) {
    snr <- snr_values[si]
    for (mc in seq_len(NUM_MC_SINGLE)) {
      seed <- 1000L * (si - 1L) + mc
      set.seed(seed + 500000L)
      sup <- sample.int(P, support_size)
      dat <- dgp_iid_snr(N, P, sup, rep(1, support_size), snr, seed = seed)
      res <- run_biobank(dat$X, matrix(dat$y, N, 1L), seed = -1L, R_boot = 500L)
      maps <- .accumulate(maps, res$statistics[1L, ], res$selected_indices[[1L]],
                          dat$true_support, si, tFDR)
    }
    cat(sprintf("  SNR %.2f -- completed %d runs.\n", snr, NUM_MC_SINGLE))
  }
  .finalize_and_save(maps, snr_values, NUM_MC_SINGLE, NUM_MC_SINGLE,
                     "d04_biobank_screen_trex_mc_snr_n300_p1000_s10_inmem")
}


# ==============================================================================
# Part 4: MC SNR sweep — multiple phenotypes
# ==============================================================================

part_4 <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Part 4: Biobank Screen-TRex — MC SNR sweep, multiple phenotypes (in-memory)\n")
  cat(sprintf("n = %d, p = %d, q = 5, s = 5, num_MC = %d\n", N, P, NUM_MC_MULTI))
  cat(strrep("=", 70L), "\n", sep = "")
  q <- 5L; support_size <- 5L; tFDR <- 0.10
  snr_values <- c(0.1, 0.5, 1.0, 2.0, 5.0)
  ns <- length(snr_values)
  maps <- .new_maps(ns)

  for (si in seq_along(snr_values)) {
    snr <- snr_values[si]
    for (mc in seq_len(NUM_MC_MULTI)) {
      x_seed <- 5000L * (si - 1L) + mc
      set.seed(x_seed); X <- matrix(rnorm(N * P), N, P)   # shared X per MC run
      true_supports <- vector("list", q); Y <- matrix(0, N, q)
      for (k in seq_len(q)) {
        set.seed(77L * 100000L + q * (NUM_MC_MULTI * (si - 1L) + (mc - 1L)) + k)
        sup <- sample.int(P, support_size)
        true_supports[[k]] <- sort(sup)
        beta <- numeric(P); beta[sup] <- 1
        Y[, k] <- .make_response(X, beta, snr)
      }
      res <- run_biobank(X, Y, seed = -1L, R_boot = 500L)
      for (k in seq_len(q)) {
        maps <- .accumulate(maps, res$statistics[k, ], res$selected_indices[[k]],
                            true_supports[[k]], si, tFDR)
      }
    }
    cat(sprintf("  SNR %.2f -- completed %d runs (%d phenotypes each).\n",
                snr, NUM_MC_MULTI, q))
  }
  total <- q * NUM_MC_MULTI
  .finalize_and_save(maps, snr_values, NUM_MC_MULTI, total,
                     "d04_biobank_screen_trex_mc_multi_n300_p1000_q5_s5_inmem")
}


#' Normalize accumulators and write the table.
#' @param cond_total denominator for est.FDR / usage (total phenotype screenings).
.finalize_and_save <- function(maps, snr_values, num_MC, cond_total, stem) {
  ns <- length(snr_values)
  for (m in BIOBANK_METHODS) {
    for (i in seq_len(ns)) {
      if (maps$usage[[m]][i] > 0) {
        maps$fdp[[m]][i] <- maps$fdp[[m]][i] / maps$usage[[m]][i]   # conditional
        maps$tpp[[m]][i] <- maps$tpp[[m]][i] / maps$usage[[m]][i]
      }
      maps$est[[m]][i]   <- maps$est[[m]][i]   / cond_total          # unconditional
      maps$usage[[m]][i] <- maps$usage[[m]][i] / cond_total          # fraction
    }
  }
  save_and_print_biobank_mc(num_MC, OUT_DIR, stem, snr_values, BIOBANK_METHODS,
                            maps$fdp, maps$tpp, maps$est, maps$usage)
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_1) part_1()
if (run_part_2) part_2()
if (run_part_3) part_3()
if (run_part_4) part_4()
cat("\ndemo_trex_scr_04_biobank_inmem.R complete.\n")
