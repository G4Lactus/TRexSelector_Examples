# ==============================================================================
# demo_trex_da_08_groups.R
# ==============================================================================
#
# DA-TRex MC simulation for the prior-groups DGP. The predictors carry a known
# 3-level hierarchical group structure (groups of 10 / 50 / 250); the prior-
# groups DA method uses that structure directly. Data is drawn from a multi-level
# latent-factor model with Toeplitz leaf blocks.
#
#  Part 1: SNR sweep (n=300, p=1000, s=10, Random support).
#
# The method is selected by supplying a non-empty `prior_groups` to
# `trex_da_control()` (the core routes on prior_groups being non-empty; the
# `da_method` value is then irrelevant). Mirrors
# cpp/.../demo_trex_da_08_mc_sim_groups and the Python port.
# ==============================================================================

library(TRexSelectorNeo)
library(parallel)

num_cores <- 6L

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)

source(file.path(this_dir_, "..", "..", "support_generators.R"))   # make_support_random
source(file.path(this_dir_, "..", "trex_da_dgps.R"))
source(file.path(this_dir_, "..", "..", "simulation_utils.R"))     # .print_table


# ==============================================================================
# Base parameters
# ==============================================================================
BASE <- list(
  n         = 300L,
  p         = 1000L,
  s         = 10L,
  amplitude = 1.0,
  tFDR      = 0.2,
  K         = 20L,
  num_MC    = 200L,
  seed      = 2026L,
  max_gap   = 20L
)

# 3-level hierarchy: groups of 10, 50, 250 (fine -> coarse), 0-based ids.
group_structure <- function(p) {
  list((seq_len(p) - 1L) %/% 10L,
       (seq_len(p) - 1L) %/% 50L,
       (seq_len(p) - 1L) %/% 250L)
}
RHO_LEVELS <- c(0.55, 0.25, 0.10)
PHI_LEAF   <- 0.50


# ==============================================================================
# MC runner (data seed = base + trial, selector seed = -1)
# ==============================================================================
.run_mc_groups <- function(snr) {
  groups <- group_structure(BASE$p)

  one_trial <- function(mc) {
    trial_seed <- BASE$seed + mc - 1L
    support    <- make_support_random(BASE$s, BASE$p, trial_seed)
    dat <- dgp_groups_toeplitz_leaf(
      n = BASE$n, p = BASE$p, support = support, amplitude = BASE$amplitude,
      groups = groups, rho_levels = RHO_LEVELS, phi = PHI_LEAF,
      snr = snr, seed = trial_seed)

    sel <- TRexDASelector$new(
      dat$X, dat$y, tFDR = BASE$tFDR, seed = -1L, verbose = FALSE,
      da_control = trex_da_control(prior_groups = groups, rho_grid_labels = RHO_LEVELS),
      control    = trex_control(solver = "TLARS", K = BASE$K))
    sel$select()

    ts <- which(dat$beta != 0)
    c(FDP = compute_fdp(sel$selected_indices, ts),
      TPP = compute_tpp(sel$selected_indices, ts))
  }

  res <- if (.Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(BASE$num_MC), one_trial, mc.cores = num_cores)
  } else {
    lapply(seq_len(BASE$num_MC), one_trial)
  }
  M <- do.call(rbind, res)
  list(mean_FDP = mean(M[, "FDP"]), mean_TPP = mean(M[, "TPP"]),
       sd_FDP   = stats::sd(M[, "FDP"]), sd_TPP = stats::sd(M[, "TPP"]))
}


# ==============================================================================
# Part 1: SNR sweep
# ==============================================================================
if (TRUE) {
  snr_grid <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0)

  cat(strrep("=", 70), "\n")
  cat(sprintf(
    "  Prior-groups + Toeplitz leaf -- SNR sweep  |  n=%d  p=%d  s=%d  %d MC\n",
    BASE$n, BASE$p, BASE$s, BASE$num_MC))
  cat(sprintf("  3-level hierarchy (groups of 10 / 50 / 250), rho_levels={%s}, phi_leaf=%.2f\n",
              paste(RHO_LEVELS, collapse = ", "), PHI_LEAF))
  cat(strrep("=", 70), "\n\n")

  snr_results <- lapply(snr_grid, function(snr_val) {
    cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, BASE$num_MC))
    c(list(SNR = snr_val), .run_mc_groups(snr_val))
  })

  cat("\n  Results -- SNR sweep:\n")
  .print_table(snr_results, "SNR", "%-8.2f")
}

cat("\ndemo_trex_da_08_groups.R complete.\n")
