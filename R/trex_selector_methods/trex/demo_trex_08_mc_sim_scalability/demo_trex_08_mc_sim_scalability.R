# ==============================================================================
# demo_trex_08_mc_sim_scalability.R
# ==============================================================================
#
# Part 7 of the classical T-Rex selector demos. This part demonstrates
# a Monte Carlo (MC) simulation to benchmark the scalability of the
# T-Rex selector's implementation.
#
# We compare the in-memory execution against the chunked, memory-mapped
# out-of-core execution over an exponentially increasing grid of n and p.
#
# Demo content:
# ----------------------------------------------------------------
#
# - Defines a grid of (n, p) scenarios targeting specific memory footprints
#   for X.
# - For each scenario, runs the classical T-Rex Selector in-memory and with
#   memory-mapped out-of-core execution.
# - Measures execution time and peak memory usage for both approaches.
# - Results are saved to a CSV file and printed to the console.
#
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
source(file.path(this_dir_, "..", "trex_sim_utils.R"))

# =============================================================================
# Global Simulation Parameters
# =============================================================================
OUT_DIR <- file.path(this_dir_, "simulation_results", "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
run_part_7 <- TRUE


# ==============================================================================
# Benchmark Configuration
# ==============================================================================

# Define exponential (n, p) grid targeting specific memory footprints for X
# (assuming 8 bytes per double precision float)
SCENARIOS <- list(
  list(n = 5000L,  p = 25000L,  label = "1 GB"),
  list(n = 15000L, p = 75000L,  label = "9 GB"),
  list(n = 30000L, p = 150000L, label = "36 GB"),
  list(n = 50000L, p = 250000L, label = "100 GB"),
  list(n = 80000L, p = 400000L, label = "256 GB")
)

S_TRUE   <- 10L
SNR      <- 1.0
SEED     <- 42L

# ==============================================================================
# Execution Helpers
# ==============================================================================

run_in_memory <- function(n, p, s, snr, seed, ctrl) {
  gc(reset = TRUE)
  start_time <- Sys.time()

  cat("    -> Generating in-memory data...\n")
  dat <- dgp_gauss_snr(n, p, seq_len(s), snr = snr, seed = seed)

  cat("    -> Running T-Rex...\n")
  sel <- TRexSelector$new(X = dat$X, y = dat$y, tFDR = 0.1, control = ctrl)
  sel$select()

  end_time <- Sys.time()
  mem_mb   <- sum(gc()[, 6L]) # Column 6 is "Max used (MB)" in R

  list(
    time_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
    mem_mb   = mem_mb,
    status   = "OK"
  )
}


run_mmap <- function(n, p, s, snr, seed, ctrl) {
  gc(reset = TRUE)
  start_time <- Sys.time()

  x_path <- tempfile(fileext = ".dat")
  on.exit(unlink(x_path), add = TRUE)

  cat(sprintf("    -> Generating chunked out-of-core data to %s...\n", x_path))
  dat <- dgp_chunked_gauss_snr(n, p, seq_len(s), x_path,
                               chunk_cols = 1000L, snr = snr, seed = seed)

  cat("    -> Running memory-mapped T-Rex...\n")
  sel <- TRexSelector$new(X = dat$X, y = dat$y, tFDR = 0.1, control = ctrl)
  sel$select()

  end_time <- Sys.time()
  mem_mb   <- sum(gc()[, 6L])

  list(
    time_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
    mem_mb   = mem_mb,
    status   = "OK"
  )
}


# ==============================================================================
# Main Benchmark Loop
# ==============================================================================

run_scalability_benchmark <- function() {
  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Scalability Benchmark: In-Memory vs. Memory-Mapped\n")
  cat(strrep("=", 70L), "\n\n")

  # We'll use the TLARS solver for the benchmark
  ctrl_inmem <- do.call(trex_control, list(solver = "TLARS"))
  ctrl_mmap  <- make_mmap_trex_ctrl(solver = "TLARS")

  results <- list()

  for (i in seq_along(SCENARIOS)) {
    scen <- SCENARIOS[[i]]
    n <- scen$n
    p <- scen$p
    lbl <- scen$label

    cat(sprintf("Scenario %d: n = %d, p = %d (~%s of raw data)\n", i, n, p, lbl))
    cat(strrep("-", 50L), "\n")

    # 1. In-Memory
    cat("  [In-Memory]\n")
    res_inmem <- tryCatch({
      run_in_memory(n, p, S_TRUE, SNR, SEED, ctrl_inmem)
    }, error = function(e) {
      cat("    -> Failed! Error:", conditionMessage(e), "\n")
      list(time_sec = NA_real_, mem_mb = NA_real_, status = "OOM")
    })
    cat(sprintf("    Status: %s | Time: %.1f s | Max RAM: %.0f MB\n",
                res_inmem$status, res_inmem$time_sec, res_inmem$mem_mb))

    # 2. Memory-Mapped
    cat("  [Memory-Mapped]\n")
    res_mmap <- tryCatch({
      run_mmap(n, p, S_TRUE, SNR, SEED, ctrl_mmap)
    }, error = function(e) {
      cat("    -> Failed! Error:", conditionMessage(e), "\n")
      list(time_sec = NA_real_, mem_mb = NA_real_, status = "ERROR")
    })
    cat(sprintf("    Status: %s | Time: %.1f s | Max RAM: %.0f MB\n",
                res_mmap$status, res_mmap$time_sec, res_mmap$mem_mb))

    cat("\n")

    # Append to results
    results <- c(results, list(
      data.frame(
        scenario   = lbl,
        n          = n,
        p          = p,
        type       = "In-Memory",
        time_sec   = res_inmem$time_sec,
        mem_mb     = res_inmem$mem_mb,
        status     = res_inmem$status,
        stringsAsFactors = FALSE
      ),
      data.frame(
        scenario   = lbl,
        n          = n,
        p          = p,
        type       = "Memory-Mapped",
        time_sec   = res_mmap$time_sec,
        mem_mb     = res_mmap$mem_mb,
        status     = res_mmap$status,
        stringsAsFactors = FALSE
      )
    ))
  }

  df_res <- do.call(rbind, results)

  csv_path <- file.path(OUT_DIR, "d04_scalability_benchmark.csv")
  utils::write.csv(df_res, csv_path, row.names = FALSE, quote = FALSE)
  cat(sprintf("\n[Info] Benchmark finished. Results saved to: %s\n\n", csv_path))

  print(df_res)
}


# ==============================================================================
# Main
# ----------
# Run part 7: Scalability Benchmark
# ==============================================================================
if (run_part_7) run_scalability_benchmark()
cat("\ndemo_trex_08_mc_sim_scalability.R complete.\n")
