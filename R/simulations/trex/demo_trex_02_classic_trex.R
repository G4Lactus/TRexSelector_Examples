# ==============================================================================
# demo_trex_02_classic_trex.R
# ==============================================================================
#
# Mirror of cpp/trex_selector_methods/trex/demo_trex_02_classic_trex_mmap.cpp
#
# The C++ file demonstrates two mmap usage patterns. This R file mirrors the
# same statistical scenarios; the memory-mapping aspect is handled as follows:
#
#  Part A: In-memory X + use_memory_mapping = TRUE.
#            Enables internal D-mmap and solver serialization within
#            TRexSelector (via trex_control(use_memory_mapping = TRUE)).
#            Exact mirror of C++ Demo A.
#
#  Part B: Fully memory-mapped X + use_memory_mapping = TRUE.
#            In C++ this uses SyntheticDataMapped (X never fully in RAM).
#            The R package TRexSelector supports MemoryMappedMatrix objects
#            for X (via convert_to_memory_mapped / mmap_matrix), but creating
#            a memory-mapped version of a generated matrix requires writing it
#            to disk first.  This part mirrors the statistical scenario
#            (same n, p, support, snr, ctrl) while demonstrating the R mmap
#            interface for completeness.
#
# ==============================================================================

library(TRexSelector)

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) "."
)
source(file.path(this_dir_, "dgp_gauss_snr.R"))
source(file.path(this_dir_, "..", "simulation_utils.R"))


# ==============================================================================
# Parameters
# ==============================================================================

OUT_DIR <- file.path(this_dir_, "simulation_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Run flags
run_part_a <- TRUE
run_part_b <- TRUE

# Shared demo config (mirrors make_mmap_demo_config() in C++)
# C++ true_support (0-based): {27, 149, 398, 420, 4}  →  R (1-based):
DEMO_CFG <- list(
  true_support = c(28L, 150L, 399L, 421L, 5L),
  true_coefs   = c(1, 1, 1, 1, 1),   # rnd_coef = FALSE
  snr          = 1.0,
  tFDR         = 0.1
)

# Shared TRex control (mirrors make_mmap_trex_ctrl() in C++)
# use_memory_mapping = TRUE enables D-mmap + solver serialization.
MMAP_CTRL_BASE <- list(
  K                        = 20L,
  max_dummy_multiplier     = 10L,
  use_max_T_stop           = TRUE,
  lloop_strategy           = "HCONCAT",
  tloop_stagnation_stop    = TRUE,
  use_memory_mapping       = TRUE
  # method = "TLARS" added below
)


# ==============================================================================
# Shared helper
# ==============================================================================

#' Print and save single-run results.
.print_and_save_d02 <- function(stem, dat, selector, tFDR) {
  true_support <- dat$true_support
  selected     <- selector$selected_indices
  n_sel        <- length(selected)
  n_tp         <- length(intersect(selected, true_support))
  n_fp         <- n_sel - n_tp
  tpp          <- TRexSelector::compute_tpp(selected, true_support)
  fdp          <- TRexSelector::compute_fdp(selected, true_support)

  cat(strrep("-", 70L), "\n")
  cat(sprintf("  Selected indices: {%s}\n", paste(selected, collapse = ", ")))
  cat(sprintf("  FDP = %.4f  |  TPP = %.4f  (target tFDR <= %.2f)\n",
              fdp, tpp, tFDR))
  cat(sprintf("  T_stop = %d,  L = %d\n", selector$T_stop, selector$L))
  cat(strrep("-", 70L), "\n\n")

  # Save per-variable CSV (mirrors save_selection_csv() in C++)
  p          <- ncol(dat$X)
  phi_prime  <- selector$phi_prime
  is_sel     <- seq_len(p) %in% selected
  is_true    <- seq_len(p) %in% true_support
  csv_df     <- data.frame(
    variable_index = seq_len(p) - 1L,   # 0-based to match C++ output
    phi_prime      = round(phi_prime, 6L),
    selected       = as.integer(is_sel),
    true_positive  = as.integer(is_true)
  )
  csv_path <- file.path(OUT_DIR, paste0(stem, ".csv"))
  utils::write.csv(csv_df, csv_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("[Info] CSV saved to: %s\n\n", csv_path))
}


# ==============================================================================
# Part A: In-memory X + use_memory_mapping = TRUE
# ==============================================================================
# Mirrors demo_TRexSelector_d_mmap_solver_serial() in C++.

trx_02_part_a <- function() {

  run_single <- function(high_dim) {
    cfg       <- DEMO_CFG
    cfg$n     <- if (high_dim) 1000L else 5000L
    cfg$p     <- if (high_dim) 5000L else 1000L
    dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"

    cat("\n", strrep("=", 70L), "\n", sep = "")
    cat("Part A: In-Memory X  +  use_memory_mapping = TRUE\n")
    cat("  Activates D-mmap and solver serialization inside TRexSelector.\n")
    cat(sprintf("  %s  (n = %d, p = %d)\n\n", dim_label, cfg$n, cfg$p))

    ctrl <- do.call(trex_control, c(list(method = "TLARS"), MMAP_CTRL_BASE))

    cat("Generating synthetic data (in-memory) ...\n")
    dat <- dgp_gauss_snr(cfg$n, cfg$p, cfg$true_support, coefs = cfg$true_coefs,
                         snr = cfg$snr, seed = 58L)

    cat("Running T-Rex Selector ...\n\n")
    selector <- TRexSelector$new(
      X       = dat$X,
      y       = dat$y,
      tFDR    = cfg$tFDR,
      seed    = 58L,
      verbose = TRUE,
      control = ctrl
    )
    selector$select()

    stem <- sprintf("d02_trex_mmap_demo_a_n%d_p%d", cfg$n, cfg$p)
    .print_and_save_d02(stem, dat, selector, cfg$tFDR)
  }

  run_single(high_dim = FALSE)
  run_single(high_dim = TRUE)
}


# ==============================================================================
# Part B: Memory-mapped X + use_memory_mapping = TRUE
# ==============================================================================
# In C++: uses SyntheticDataMapped so X is never fully in RAM.
# In R:   demonstrates the TRexSelector MemoryMappedMatrix interface by
#         writing X to a temporary mmap file before passing it to the selector.
#         use_memory_mapping = TRUE adds D-mmap + solver serialization on top.

trx_02_part_b <- function() {

  run_single <- function(high_dim) {
    cfg       <- DEMO_CFG
    cfg$n     <- if (high_dim) 1000L else 5000L
    cfg$p     <- if (high_dim) 5000L else 1000L
    dim_label <- if (high_dim) "High-dimensional (p > n)" else "Low-dimensional (n > p)"

    cat("\n", strrep("=", 70L), "\n", sep = "")
    cat("Part B: Memory-Mapped X  +  use_memory_mapping = TRUE\n")
    cat("  Full mmap pipeline: X mmap + D mmap + solver serialization.\n")
    cat(sprintf("  %s  (n = %d, p = %d)\n\n", dim_label, cfg$n, cfg$p))

    ctrl <- do.call(trex_control, c(list(method = "TLARS"), MMAP_CTRL_BASE))

    # Generate data in-memory, then convert X to a memory-mapped file.
    cat("Generating synthetic data ...\n")
    dat <- dgp_gauss_snr(cfg$n, cfg$p, cfg$true_support, coefs = cfg$true_coefs,
                         snr = cfg$snr, seed = 58L)

    # Write X to a temporary mmap backing file (RAII-style: deleted on exit)
    x_path <- tempfile(fileext = ".dat")
    on.exit(unlink(x_path), add = TRUE)

    cat(sprintf("Converting X to memory-mapped file: %s ...\n", x_path))
    X_mmap <- TRexSelector::convert_to_memory_mapped(dat$X, x_path)

    cat("Running T-Rex Selector (X mmap + D mmap + solver serial) ...\n\n")
    selector <- TRexSelector$new(
      X       = X_mmap,
      y       = dat$y,
      tFDR    = cfg$tFDR,
      seed    = 58L,
      verbose = TRUE,
      control = ctrl
    )
    selector$select()

    stem <- sprintf("d02_trex_mmap_demo_b_n%d_p%d", cfg$n, cfg$p)
    .print_and_save_d02(stem, dat, selector, cfg$tFDR)
  }

  run_single(high_dim = FALSE)
  run_single(high_dim = TRUE)
}


# ==============================================================================
# Main
# ==============================================================================

if (run_part_a) trx_02_part_a()
if (run_part_b) trx_02_part_b()

cat("\ndemo_trex_02_classic_trex.R complete.\n")
# ==============================================================================
