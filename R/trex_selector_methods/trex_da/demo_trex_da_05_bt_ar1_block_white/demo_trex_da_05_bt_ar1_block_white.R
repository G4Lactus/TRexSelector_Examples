# ==============================================================================
# demo_trex_da_06_ar1pluswhite_block_sweeps.R
# ==============================================================================
#
# Scenario 2 — AR(1)-block + white-noise DGP.
#
# DGP: dgp_ar1_block_white_snr.
#
# The total number of predictors p_total is fixed throughout all sweeps.
# Of these, p_ar = M * Q columns follow M times AR(1) blocks of size
# Q (the structured part), and p_white = p_total - M * Q columns are i.i.d. N(0,1)
# null variables (the white-noise augmentation).
# Active predictors — one per AR block — lie entirely within the AR part.
#
# This model allows p_total to remain constant while varying the block structure
# parameters M and Q, something that is impossible with the pure AR(1)-block
# DGP of Scenario 1  (demo_trex_da_05_ar1_block_sweeps.R).
#
# Each part fixes all but one structural parameter and sweeps it, with an
# outer loop over hc_dist in {"single", "complete", "average"}.
#
#  Part 1: SNR sweep    SNR in {0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}
#                       p_ar = 25, p_white = 975 fixed.
#
#  Part 2: rho sweep    rho in {0.0, 0.1, ..., 0.9}
#                       p_ar = 25, p_white = 975 fixed.
#
#  Part 3: Q sweep      Q   in {5, 10, 15, 20, 25, 30, 35, 40, 45, 50}
#                       p_ar = M*Q in {25..250}; p_white = p_total - p_ar.
#
#  Part 4: M sweep      M   in {1, 3, 5, 10, 15, 20, 30}
#                       p_ar = M*Q in {5..150}; p_white = p_total - p_ar; s = M.
#
# ==============================================================================

this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)
source(file.path(this_dir_, "..", "support_generators.R"))   # make_support_bt_one_per_block
source(file.path(this_dir_, "dgp_bt_snr.R"))                 # dgp_ar1_block_white_snr
source(file.path(this_dir_, "..", "simulation_utils.R"))

library(TRexSelectorNeo)
library(parallel)
num_cores <- 6L

# ==============================================================================
# Base parameters
# ==============================================================================
BASE <- list(
  n         = 300L,   # observations
  p_total   = 1000L,  # total predictors (fixed throughout all sweeps)
  M         = 5L,     # number of AR(1) blocks (= number of active predictors)
  Q         = 5L,     # block size; p_ar = M*Q = 25 at base setting
  amplitude = 1.0,    # signal coefficient (calibrated out by SNR)
  rho       = 0.7,    # AR(1) within-block correlation
  snr       = 2.0,    # signal-to-noise ratio
  tFDR      = 0.1,    # target FDR
  K         = 20L,    # random experiments per selector run
  num_MC    = 200L,   # MC replications
  seed      = 2026L   # base RNG seed
)

linkage_grid <- c("single", "complete", "average")


# ==============================================================================
#  Part 1: SNR sweep
#  Fixed: n, p = 1000, M = 5, Q = 5 (p_ar = 25, p_white = 975, s = 5), rho = 0.7.
#  Swept: SNR in {0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_06_bt_ar1_pl_wh_blk_01 <- function() {
    snr_grid   <- c(0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
    p_ar_p0    <- BASE$M * BASE$Q                # 25
    p_white_p0 <- BASE$p_total - p_ar_p0         # 975
    support_p0 <- make_support_bt_one_per_block(p = p_ar_p0, n_blocks = BASE$M)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)+white -- SNR sweep  |  n=%d  p_total=%d  M=%d  Q=%d  p_ar=%d  p_white=%d  s=%d  rho=%.1f  %d MC\n",
      BASE$n, BASE$p_total, BASE$M, BASE$Q, p_ar_p0, p_white_p0, BASE$M, BASE$rho, BASE$num_MC
    ))
    cat(sprintf("  SNR in {%s}\n", paste(snr_grid, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      snr_results <- lapply(snr_grid, function(snr_val) {
        cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, BASE$num_MC))
        r <- .run_mc_ar1_block_white(
          n         = BASE$n, p_ar = p_ar_p0, p_white = p_white_p0,
          support   = support_p0, amplitude = BASE$amplitude,
          n_blocks  = BASE$M, rho = BASE$rho, hc_dist = lnk,
          snr       = snr_val, tFDR = BASE$tFDR, K = BASE$K,
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
  trx_da_06_bt_ar1_pl_wh_blk_01()
}  # end Part 1
# ==============================================================================


# ==============================================================================
#  Part 2: rho sweep
#  Fixed: n, p=1000, M=5, Q=5 (p_ar=25, p_white=975, s=5), SNR=2.
#  Swept: rho in {0.0, 0.1, 0.2, ..., 0.9}.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_06_bt_ar1_pl_wh_blk_02 <- function() {
    rho_grid   <- seq(0, 0.9, by = 0.1)
    p_ar_p1    <- BASE$M * BASE$Q                # 25
    p_white_p1 <- BASE$p_total - p_ar_p1         # 975
    support_p1 <- make_support_bt_one_per_block(p = p_ar_p1, n_blocks = BASE$M)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)+white -- rho sweep  |  n=%d  p_total=%d  M=%d  Q=%d  p_ar=%d  p_white=%d  s=%d  SNR=%.1f  %d MC\n",
      BASE$n, BASE$p_total, BASE$M, BASE$Q, p_ar_p1, p_white_p1, BASE$M, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  rho in {%s}\n", paste(rho_grid, collapse = ", ")))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      rho_results <- lapply(rho_grid, function(rho_val) {
        cat(sprintf("  [rho = %.1f] running %d MC trials ...\n", rho_val, BASE$num_MC))
        r <- .run_mc_ar1_block_white(
          n         = BASE$n, p_ar = p_ar_p1, p_white = p_white_p1,
          support   = support_p1, amplitude = BASE$amplitude,
          n_blocks  = BASE$M, rho = rho_val, hc_dist = lnk,
          snr       = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(rho = rho_val), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table(rho_results, "rho", "%-6.1f")
    }
  }

  # Run demo for Part 2
  trx_da_06_bt_ar1_pl_wh_blk_02()
}  # end Part 2
# ==============================================================================


# ==============================================================================
#  Part 3: Q sweep (block size / group size)
#  Fixed: n, p=1000, M=5, rho=0.7, SNR=2.
#  Swept: Q in {5, 10, 15, 20, 25, 30, 35, 40, 45, 50}.
#         p_ar = M*Q in {25..250}; p_white = 1000 - p_ar; s = M = 5 fixed.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_06_bt_ar1_pl_wh_blk_03 <- function() {
    q_grid <- seq(5L, 50L, by = 5L)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)+white -- Q sweep  |  n=%d  p_total=%d  M=%d  s=%d  rho=%.1f  SNR=%.1f  %d MC\n",
      BASE$n, BASE$p_total, BASE$M, BASE$M, BASE$rho, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  Q in {%s}  (p_ar = M*Q varies from %d to %d; p_white = %d - p_ar)\n",
                paste(q_grid, collapse = ", "),
                BASE$M * min(q_grid), BASE$M * max(q_grid), BASE$p_total))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      q_results <- lapply(q_grid, function(q_val) {
        p_ar    <- BASE$M * q_val
        p_white <- BASE$p_total - p_ar
        support <- make_support_bt_one_per_block(p = p_ar, n_blocks = BASE$M)

        cat(sprintf("  [Q = %d, p_ar = %d, p_white = %d] running %d MC trials ...\n",
                    q_val, p_ar, p_white, BASE$num_MC))
        r <- .run_mc_ar1_block_white(
          n         = BASE$n, p_ar = p_ar, p_white = p_white,
          support   = support, amplitude = BASE$amplitude,
          n_blocks  = BASE$M, rho = BASE$rho, hc_dist = lnk,
          snr       = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(Q = q_val, p_ar = p_ar, p_white = p_white), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table_multi(q_results, c("Q", "p_ar", "p_white"),
                         c(Q = "%-4d", p_ar = "%-6d", p_white = "%-8d"))
    }
  }

  # Run demo for Part 3
  trx_da_06_bt_ar1_pl_wh_blk_03()
}  # end Part 3
# ==============================================================================


# ==============================================================================
#  Part 4: M sweep (number of groups / blocks)
#  Fixed: n, p=1000, Q=5, rho=0.7, SNR=2.
#  Swept: M in {1, 3, 5, 10, 15, 20, 30}.
#         p_ar = M*Q in {5..150}; p_white = 1000 - p_ar; s = M varies.
#  Outer loop: hc_dist in {"single", "complete", "average"}.
# ==============================================================================
if (TRUE) {
  trx_da_06_bt_ar1_pl_wh_blk_04 <- function() {
    m_grid <- c(1L, 3L, 5L, 10L, 15L, 20L, 30L)

    cat(strrep("=", 70), "\n")
    cat(sprintf(
      "  AR(1)+white -- M sweep  |  n=%d  p_total=%d  Q=%d  rho=%.1f  SNR=%.1f  %d MC\n",
      BASE$n, BASE$p_total, BASE$Q, BASE$rho, BASE$snr, BASE$num_MC
    ))
    cat(sprintf("  M in {%s}  (p_ar = M*%d varies from %d to %d; p_white = %d - p_ar; s = M)\n",
                paste(m_grid, collapse = ", "), BASE$Q,
                min(m_grid) * BASE$Q, max(m_grid) * BASE$Q, BASE$p_total))
    cat(strrep("=", 70), "\n\n")

    for (lnk in linkage_grid) {
      cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))

      m_results <- lapply(m_grid, function(m_val) {
        p_ar    <- m_val * BASE$Q
        p_white <- BASE$p_total - p_ar
        support <- make_support_bt_one_per_block(p = p_ar, n_blocks = m_val)

        cat(sprintf("  [M = %d, s = %d, p_ar = %d, p_white = %d] running %d MC trials ...\n",
                    m_val, m_val, p_ar, p_white, BASE$num_MC))
        r <- .run_mc_ar1_block_white(
          n         = BASE$n, p_ar = p_ar, p_white = p_white,
          support   = support, amplitude = BASE$amplitude,
          n_blocks  = m_val, rho = BASE$rho, hc_dist = lnk,
          snr       = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC    = BASE$num_MC, seed = BASE$seed,
          n_cores   = num_cores
        )
        c(list(M = m_val, s = m_val, p_ar = p_ar, p_white = p_white), r)
      })

      cat(sprintf("\n  hc_dist = \"%s\"\n", lnk))
      .print_table_multi(m_results, c("M", "s", "p_ar", "p_white"),
                         c(M = "%-4d", s = "%-4d", p_ar = "%-6d", p_white = "%-8d"))
    }
  }

  # Run demo for Part 4
  trx_da_06_bt_ar1_pl_wh_blk_04()
}  # end Part 4
# ==============================================================================

cat("\nAR(1)-block + white-noise parameter sweeps (Scenario 2) complete.\n")
