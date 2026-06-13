# ==============================================================================
# demo_trex_01_gauss.R
# ==============================================================================
#
# Classical T-Rex Selector demo on the package-shipped Gauss dataset.
#
#  Part 1: Single-run demo.
#            - Loads the built-in Gauss_data (n=50, p=100, s=3).
#            - Runs the classical T-Rex selector (TLARS solver).
#            - Reports calibration details, selected support, TPP and FDP.
#
# ==============================================================================

library(TRexSelector)


# ==============================================================================
# Parameters
# ==============================================================================

PARAMS <- list(
  tFDR = 0.2,   # target FDR level
  K    = 20L,   # number of random experiments
  seed = 1234L  # RNG seed
)

run_part_1 <- TRUE


# ==============================================================================
# Local helpers
# ==============================================================================

#' Print a formatted single-run result block.
.print_result <- function(scenario_name, dat, selector, tFDR) {

  true_support <- which(dat$support)
  s            <- length(true_support)
  p            <- ncol(dat$X)
  n_sel        <- length(selector$selected_indices)
  n_tp         <- length(intersect(selector$selected_indices, true_support))
  n_fp         <- n_sel - n_tp
  tpp          <- compute_tpp(selector$selected_indices, true_support)
  fdp          <- compute_fdp(selector$selected_indices, true_support)

  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s\n", scenario_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Data:  n = %d,  p = %d,  s = %d,  tFDR = %.2f\n",
              nrow(dat$X), p, s, tFDR))
  cat(sprintf("         True support (1-based):   {%s}\n",
              paste(true_support, collapse = ", ")))
  cat(sprintf("         Active coefficients:       {%s}\n",
              paste(dat$beta[true_support], collapse = ", ")))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Calibration:  T_stop = %d,  num_dummies = %d  (L = %d)\n",
              selector$T_stop, selector$L * p, selector$L))
  cat(sprintf("  Selection:    %d selected  |  TP = %d  FP = %d\n",
              n_sel, n_tp, n_fp))
  cat(sprintf("  Rates:        TPP = %.3f  |  FDP = %.3f  (target tFDR <= %.2f)\n",
              tpp, fdp, tFDR))
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}


# ==============================================================================
# Part 1: Single-run demo
# ==============================================================================

if (run_part_1) {
  trx_01_gauss_01 <- function() {

    # Header
    # --------------------------------------------------
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("Gauss Demo  |  Part 1: Single-run\n")
    cat(sprintf("tFDR = %.2f,  K = %d,  seed = %d\n",
                PARAMS$tFDR, PARAMS$K, PARAMS$seed))
    cat(strrep("=", 70), "\n\n")
    # --------------------------------------------------

    # Load built-in dataset
    data("Gauss_data")
    X <- Gauss_data$X
    y <- c(Gauss_data$y)

    cat(sprintf("[Part 1] Gauss_data:  n = %d,  p = %d\n", nrow(X), ncol(X)))
    cat(sprintf("[Part 1] True support (1-based): {%s},  beta_active = %.1f\n\n",
                paste(which(Gauss_data$support), collapse = ", "),
                unique(Gauss_data$beta[Gauss_data$support])))

    # Instantiate and run classical T-Rex selector (TLARS solver)
    cat("[Part 1] Running classical T-Rex selector (TLARS) ...\n\n")
    ctrl     <- trex_control(method = "TLARS", K = PARAMS$K)
    selector <- TRexSelector$new(
      X       = X,
      y       = y,
      tFDR    = PARAMS$tFDR,
      seed    = PARAMS$seed,
      verbose = TRUE,
      control = ctrl
    )
    selector$select()

    .print_result(
      "Part 1 \u2014 Gauss data  [classical T-Rex, TLARS]",
      Gauss_data, selector, PARAMS$tFDR
    )
  }

  trx_01_gauss_01()

}  # end Part 1
# ==============================================================================

cat("\nGauss demo complete.\n")
