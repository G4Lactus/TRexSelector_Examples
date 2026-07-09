# ==============================================================================
# demo_ts_09_tmp.R
# ==============================================================================
#
# Demonstration of the T-MP (Terminating Matching Pursuit) solver.
#
# Mirrors cpp/tsolvers/demo_ts_09_tmp/demo_ts_09_tmp.cpp.
#
# Parts:
#   Demo 1: Basic T-MP with early stopping (internal normalization),
#           high-dimensional (p > n) and low-dimensional (n > p) regimes.
#   Demo 2: External L2 normalization (LpNormScaler) + centered y,
#           solver run with normalize = FALSE, intercept = FALSE.
#   Demo 3: Serialization & warm start — save a partial path, reload,
#           continue, and compare against a reference run.
#   Demo 4: Memory-mapped X via convert_to_memory_mapped(). The R binding
#           accepts an mmap X (D and y stay in memory); cpp maps all three
#           and generates p = 500,000 on disk — scaled to p = 5,000 here.
#
# ==============================================================================

library(TRexSelectorNeo)

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

source(file.path(this_dir_, "..", "ts_demo_utils.R"))


# ==============================================================================
# Demo 1: Basic T-MP with Early Stopping
# ==============================================================================

demo_early_stopping <- function(high_dim, rnd_coef, T_stop) {
  print_section_header("Demo 1: Basic T-MP with Early Stopping")

  n <- if (high_dim) 1000L else 5000L
  p <- if (high_dim) 5000L else 1000L
  num_dummies <- 10L * p

  # cpp true_support is 0-based {27, 149, 398, 420, 4} -> 1-based here
  true_support <- c(28L, 150L, 399L, 421L, 5L)
  true_coefs <- if (rnd_coef) c(-0.4, -0.2, -0.8, 1.1, 2.5) else rep(1, 5)
  snr <- 1.0

  cat(if (high_dim) "High-dimensional (p > n)\n" else "Low-dimensional (n > p)\n")
  print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

  dat <- gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                            snr, seed = 42)

  solver <- TMP_Solver$new(dat$X, dat$D, dat$y, normalize = TRUE,
                           intercept = TRUE, verbose = TRUE)
  solver$execute_step(T_stop, early_stop = TRUE)

  print_selection(solver, true_support, p)
  print_selection_quality(solver, true_support, p)
  cat("\n\n")
}


# ==============================================================================
# Demo 2: External Normalization (LpNormScaler)
# ==============================================================================

demo_with_external_normalizer <- function(high_dim, rnd_coef, T_stop) {
  print_section_header("Demo 2: T-MP with External Normalization")

  n <- if (high_dim) 1000L else 5000L
  p <- if (high_dim) 5000L else 1000L
  num_dummies <- 10L * p

  # cpp true_support is 0-based {4, 27, 149, 398, 420} -> 1-based here
  true_support <- c(5L, 28L, 150L, 399L, 421L)
  true_coefs <- if (rnd_coef) c(2.5, -0.4, -0.2, -0.8, 1.1) else rep(1, 5)
  snr <- 1.0

  cat(if (high_dim) "High-dimensional (p > n)\n" else "Low-dimensional (n > p)\n")
  print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

  dat <- gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                            snr, seed = 42)

  # External normalization: center + unit-L2 columns for X and D, center y.
  cat("Applying external L2 normalization...\n")
  l2scaler_X <- LpNormScaler$new(norm_type = 2, center = TRUE)
  l2scaler_X$fit(dat$X)
  l2scaler_X$transform_inplace(dat$X)

  l2scaler_D <- LpNormScaler$new(norm_type = 2, center = TRUE)
  l2scaler_D$fit(dat$D)
  l2scaler_D$transform_inplace(dat$D)

  dat$y <- dat$y - mean(dat$y)

  solver <- TMP_Solver$new(dat$X, dat$D, dat$y, normalize = FALSE,
                           intercept = FALSE, verbose = TRUE)
  solver$execute_step(T_stop, early_stop = TRUE)

  print_selection(solver, true_support, p)
  print_selection_quality(solver, true_support, p)
  cat("\n\n")
}


# ==============================================================================
# Demo 3: Serialization & Warm Start
# ==============================================================================

demo_serialization <- function() {
  print_section_header("Demo 3: T-MP Serialization & Warm Start")

  n <- 100L
  p <- 50L
  num_dummies <- p
  T_stop_final <- 15L
  T_stop_partial <- 7L
  snr <- 1.0

  # cpp true_support is 0-based {10, 25, 40} -> 1-based here
  true_support <- c(11L, 26L, 41L)
  true_coefs <- c(2.5, -1.8, 3.2)

  print_demo_config(n, p, num_dummies, T_stop_final, true_support,
                    true_coefs, snr)

  dat <- gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                            snr, seed = 42)

  # Reference solver: run until T_stop_final.
  solver_ref <- TMP_Solver$new(dat$X, dat$D, dat$y, normalize = TRUE,
                               intercept = TRUE, verbose = FALSE)
  solver_ref$execute_step(T_stop_final, early_stop = TRUE)
  cat(sprintf("\u2713 Reference completed with %d steps\n\n",
              solver_ref$get_num_steps()))

  filename <- file.path(this_dir_, "tmp_checkpoint.bin")

  # STEP 1: Run partial path and save checkpoint.
  solver1 <- TMP_Solver$new(dat$X, dat$D, dat$y, normalize = TRUE,
                            intercept = TRUE, verbose = FALSE)
  solver1$execute_step(T_stop_partial, early_stop = TRUE)
  cat(sprintf("Checkpoint at partial stop: %d steps\n", solver1$get_num_steps()))
  solver1$save(filename)
  cat(sprintf("\u2713 Checkpoint saved at '%s'\n", filename))

  # STEP 2: Load from checkpoint and continue.
  # The R binding hydrates an existing instance (cpp: static
  # load(filename, X, D)); $new() re-binds the data views.
  solver2 <- TMP_Solver$new(dat$X, dat$D, dat$y, normalize = TRUE,
                            intercept = TRUE, verbose = FALSE)
  solver2$load(filename)
  cat(sprintf("Loaded from checkpoint: %d steps\n", solver2$get_num_steps()))
  solver2$execute_step(T_stop_final, early_stop = TRUE)

  cat("\nCOMPARISON:\n")
  cat(sprintf("Reference steps: %d\n", solver_ref$get_num_steps()))
  cat(sprintf("Reloaded steps:  %d\n", solver2$get_num_steps()))
  cat(sprintf("RSS diff:   %g\n",
              abs(tail(solver_ref$get_rss(), 1) - tail(solver2$get_rss(), 1))))
  cat(sprintf("R2 diff:    %g\n",
              abs(tail(solver_ref$get_r2(), 1) - tail(solver2$get_r2(), 1))))
  paths_match <- identical(unlist(solver_ref$get_actions()),
                           unlist(solver2$get_actions()))
  cat(if (paths_match) "\u2713 Paths match\n" else "\u2717 Paths differ!\n")

  file.remove(filename)
  cat("\u2713 Checkpoint file removed\n")
  cat("\n\n")
}


# ==============================================================================
# Demo 4: Memory-Mapped Data
# ==============================================================================

demo_memory_mapped <- function(high_dim, rnd_coef, T_stop) {
  print_section_header("Demo 4: T-MP with Memory-Mapped Data")

  n <- if (high_dim) 1000L else 5000L
  p <- if (high_dim) 5000L else 1000L  # cpp: 500,000 on disk; scaled here
  num_dummies <- 10L * p
  snr <- 1.0

  # cpp true_support is 0-based {4, 27, 149, 398, 420} -> 1-based here
  true_support <- c(5L, 28L, 150L, 399L, 421L)
  true_coefs <- if (rnd_coef) c(2.5, -0.4, -0.2, -0.8, 1.1) else rep(1, 5)

  print_demo_config(n, p, num_dummies, T_stop, true_support, true_coefs, snr)

  cat("Generating memory-mapped data...\n")
  X_file <- file.path(this_dir_, "demo_tmp_X.bin")

  dat <- gen_synthetic_data(n, p, num_dummies, true_support, true_coefs,
                            snr, seed = 42)

  on.exit({
    cat("\nCleaning up files...\n")
    if (file.exists(X_file)) file.remove(X_file)
    cat("\u2713 All files removed\n\n")
  }, add = TRUE)

  # Write X to disk and pass the mmap view to the solver. The R binding
  # accepts an mmap_matrix for X only; D and y remain in memory (cpp maps
  # all three via SyntheticDataMapped).
  X_mmap <- convert_to_memory_mapped(dat$X, X_file)
  dat$X <- NULL
  cat("\u2713 Data generated on disk\n\n")

  cat("Running T-MP on memory-mapped data...\n")
  solver <- TMP_Solver$new(X_mmap, dat$D, dat$y, normalize = TRUE,
                           intercept = TRUE, verbose = TRUE)
  solver$execute_step(T_stop, early_stop = TRUE)
  cat("\u2713 T-MP completed\n\n")

  print_selection(solver, true_support, p)
  print_selection_quality(solver, true_support, p)
}


# ==============================================================================
# Main
# ==============================================================================

print_section_header("T-MP Demo Suite")

# Demo 1: Basic usage with internal normalization
demo_early_stopping(high_dim = TRUE, rnd_coef = FALSE, T_stop = 10L)
demo_early_stopping(high_dim = FALSE, rnd_coef = FALSE, T_stop = 10L)

# Demo 2: External normalization
demo_with_external_normalizer(high_dim = FALSE, rnd_coef = FALSE, T_stop = 10L)
demo_with_external_normalizer(high_dim = TRUE, rnd_coef = FALSE, T_stop = 10L)

# Demo 3: Serialization
demo_serialization()

# Demo 4: Memory-mapped
demo_memory_mapped(high_dim = TRUE, rnd_coef = FALSE, T_stop = 10L)

print_section_header("\u2713 All demos completed successfully")
