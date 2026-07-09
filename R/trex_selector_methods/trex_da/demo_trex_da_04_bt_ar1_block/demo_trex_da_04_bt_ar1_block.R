# ==============================================================================
# demo_trex_da_05_ar1_block_sweeps.R
# ==============================================================================
#
# AR(1)-block DGP: block-diagonal Toeplitz covariance.
#
# DGP: dgp_ar1_block_snr.
# The M AR(1) blocks are the entire design matrix.
# p = M * Q grows naturally with the sweep parameters.
#
# Each part fixes all but one structural parameter and sweeps it, with an
# outer loop over hc_dist in {"single", "complete", "average"}.
#
#  Part 1: SNR sweep    SNR in {0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}
#
#  Part 2: rho sweep    rho in {0.0, 0.1, ..., 0.9}
#
#  Part 3: Q sweep      Q   in {5, 10, 15, 20, 25, 30, 35, 40, 45, 50}
#                       p = M * Q varies {25..250}; s = M = 5 fixed.
#
#  Part 4: M sweep      M   in {1, 3, 5, 10, 15, 20, 30}
#                       p = M * Q and s = M both vary.
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
source(file.path(this_dir_, "..", "..", "support_generators.R"))  # make_support_bt_one_per_block
source(file.path(this_dir_, "..", "trex_da_dgps.R"))
source(file.path(this_dir_, "..", "..", "simulation_utils.R"))

# ==============================================================================
# Base parameters
# ==============================================================================
BASE <- list(
  n         = 150L,   # observations
  M         = 5L,     # number of AR(1) blocks
  Q         = 5L,     # block size (group size); p = M * Q
  amplitude = 1.0,    # signal coefficient (calibrated out by SNR)
  rho       = 0.7,    # AR(1) within-block correlation parameter
  snr       = 2.0,    # signal-to-noise ratio
  tFDR      = 0.2,    # target FDR (alpha in paper)
  K         = 20L,    # random experiments per selector run
  num_MC    = 200L,   # MC replications (paper: 955)
  seed      = 2026L   # base RNG seed
)

linkage_grid <- c("single", "complete", "average")


# ==============================================================================
#  Part 1: SNR sweep
#  -------------------------------------------
#  Fixed: n, M = 5, Q = 5 (p = 25, s = 5), rho = 0.7.
#  Swept: SNR in {0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_05_bt_ar1_blswe_01 <- function() {
    # Config parameters
    snr_grid   <- c(0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
    p_p0       <- BASE$M * BASE$Q                    # p = 25
    support_p0 <- make_support_bt_one_per_block(p = p_p0, n_blocks = BASE$M)

    # Header
    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)-block -- SNR sweep  |  n=%d  M=%d  Q=%d  p=%d  s=%d  rho=%.1f  %d MC\n",
      BASE$n, BASE$M, BASE$Q, p_p0, BASE$M, BASE$rho, BASE$num_MC
    ))
    cat(sprintf("  SNR in {%s}\n", paste(snr_grid, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    # Run MC trials for each SNR value
    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      snr_results <- lapply(snr_grid, function(snr_val) {
        cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, BASE$num_MC))
        r <- .run_mc_ar1_block(
          n         = BASE$n, p = p_p0, support = support_p0,
          amplitude = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
          hc_dist   = lnk, snr = snr_val, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(SNR = snr_val), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table(snr_results, "SNR", "%-6.2f")
    }
  }

  # Run demo for Part 1
  trx_da_05_bt_ar1_blswe_01()
}  # end Part 1
# ==============================================================================



# ==============================================================================
#  Part 2: rho sweep
#  -------------------------------------------
#  Fixed: n, M = 5, Q = 5 (p = 25, s = 5), SNR = 2.
#  Swept: rho in {0.0, 0.1, 0.2, ..., 0.9}.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_05_bt_ar1_blswe_02 <- function() {
    rho_grid   <- seq(0, 0.9, by = 0.1)
    p_p1       <- BASE$M * BASE$Q                    # p = 25
    support_p1 <- make_support_bt_one_per_block(p = p_p1, n_blocks = BASE$M)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)-block -- rho sweep  |  n=%d  M=%d  Q=%d  p=%d  s=%d  SNR=%.1f  %d MC\n",
      BASE$n, BASE$M, BASE$Q, p_p1, BASE$M, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  rho in {%s}\n", paste(rho_grid, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      rho_results <- lapply(rho_grid, function(rho_val) {
        cat(sprintf("  [rho = %.1f] running %d MC trials ...\n", rho_val, BASE$num_MC))
        r <- .run_mc_ar1_block(
          n         = BASE$n, p = p_p1, support = support_p1,
          amplitude = BASE$amplitude, n_blocks = BASE$M, rho = rho_val,
          hc_dist   = lnk, snr = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(rho = rho_val), r)
      })

      # Output table
      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table(rho_results, "rho", "%-6.1f")
    }
  }

  # Run demo for Part 2
  trx_da_05_bt_ar1_blswe_02()
}  # end Part 2
# ==============================================================================



# ==============================================================================
#  Part 3: Q sweep (block size / group size)
#  -------------------------------------------
#  Fixed: n, M = 5, rho = 0.7, SNR = 2.
#  Swept: Q in {5, 10, 15, 20, 25, 30, 35, 40, 45, 50}; p = M*Q varies.
#  s = M = 5 (one active per block) throughout.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_05_bt_ar1_blswe_03 <- function() {
    q_grid <- seq(5L, 50L, by = 5L)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)-block -- Q sweep  |  n=%d  M=%d  s=%d  rho=%.1f  SNR=%.1f  %d MC\n",
      BASE$n, BASE$M, BASE$M, BASE$rho, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  Q in {%s}  (p = M*Q varies from %d to %d)\n",
                paste(q_grid, collapse = ", "),
                BASE$M * min(q_grid), BASE$M * max(q_grid)))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      q_results <- lapply(q_grid, function(q_val) {
        p_val   <- BASE$M * q_val
        support <- make_support_bt_one_per_block(p = p_val, n_blocks = BASE$M)
        cat(sprintf("  [Q = %d, p = %d] running %d MC trials ...\n",
                    q_val, p_val, BASE$num_MC))
        r <- .run_mc_ar1_block(
          n         = BASE$n, p = p_val, support = support,
          amplitude = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
          hc_dist   = lnk, snr = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(Q = q_val, p = p_val), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table_multi(q_results, c("Q", "p"), c(Q = "%-4d", p = "%-6d"))
    }
  }

  # Run demo for Part 3
  trx_da_05_bt_ar1_blswe_03()
}  # end Part 3
# ==============================================================================



# ==============================================================================
#  Part 4: M sweep (number of groups / blocks)
#  -------------------------------------------
#  Fixed: n, Q = 5, rho = 0.7, SNR = 2.
#  Swept: M in {1, 3, 5, 10, 15, 20, 30}; p = M * Q and s = M both vary.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_05_bt_ar1_blswe_04 <- function() {
    m_grid <- c(1L, 3L, 5L, 10L, 15L, 20L, 30L)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)-block -- M sweep  |  n=%d  Q=%d  rho=%.1f  SNR=%.1f  %d MC\n",
      BASE$n, BASE$Q, BASE$rho, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  M in {%s}  (p = M*%d, s = M vary)\n",
                paste(m_grid, collapse = ", "), BASE$Q))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      m_results <- lapply(m_grid, function(m_val) {
        p_val   <- m_val * BASE$Q
        support <- make_support_bt_one_per_block(p = p_val, n_blocks = m_val)
        cat(sprintf("  [M = %d, p = %d, s = %d] running %d MC trials ...\n",
                    m_val, p_val, m_val, BASE$num_MC))
        r <- .run_mc_ar1_block(
          n         = BASE$n, p = p_val, support = support,
          amplitude = BASE$amplitude, n_blocks = m_val, rho = BASE$rho,
          hc_dist   = lnk, snr = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(M = m_val, p = p_val, s = m_val), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table_multi(m_results, c("M", "p", "s"),
                         c(M = "%-4d", p = "%-6d", s = "%-4d"))
    }
  }

  # Run demo for Part 4
  trx_da_05_bt_ar1_blswe_04()
}  # end Part 4
# ==============================================================================

cat("\nAR(1)-block parameter sweeps (Scenario 1) complete.\n")
