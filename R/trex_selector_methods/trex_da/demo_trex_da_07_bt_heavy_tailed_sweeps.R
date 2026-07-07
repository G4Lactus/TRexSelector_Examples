# ==============================================================================
# demo_trex_da_07_bt_heavy_tailed_sweeps.R
# ==============================================================================
#
# Output from this script is redirected to a log file in simulation_results/
# next to this script (see below).
# To view output in the console instead, comment out the sink() lines.
#
# Testing trex+DA+BT against Heavy-Tailed Student-t Data (nu = 3).
#
# Scenario 1: Predictors X are multivariate Student-t (nu=3).
#                     Noise is Gaussian.
# Scenario 2: Predictors X are multivariate Student-t (nu=3).
#                      Noise is Student-t (nu=3).
#
# Covariance: Block-diagonal AR(1) Toeplitz (rho = 0.8 by default).
#
# Sweeps over: SNR, rho, Q (block size), M (number of blocks).
# ==============================================================================

library(TRexSelector)
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

log_file <- file.path(this_dir_, "simulation_results",
                      "demo_trex_da_07_bt_heavy_tailed_sweeps_log.txt")
sink(log_file)
on.exit(sink(), add = TRUE)

source(file.path(this_dir_, "..", "support_generators.R"))  # make_support_bt_one_per_block
source(file.path(this_dir_, "dgp_ht_snr.R"))           # dgp_ht_block_snr
source(file.path(this_dir_, "..", "simulation_utils.R"))


# ==============================================================================
# Base parameters
# ==============================================================================
BASE <- list(
  n         = 150L,    # observations
  M         = 5L,      # number of AR(1) blocks
  Q         = 5L,      # block size (group size); p = M * Q
  amplitude = 1.0,     # signal coefficient (calibrated out by SNR)
  rho       = 0.8,     # default AR(1) within-block correlation parameter
  snr       = 2.0,     # default signal-to-noise ratio
  nu        = 3.0,     # degrees of freedom for Student-t
  tFDR      = 0.2,     # target FDR
  K         = 20L,     # random experiments per selector run
  num_MC    = 200L,    # MC replications
  seed      = 2026L    # base RNG seed
)

linkage_grid  <- c("single", "complete", "average")
scenario_grid <- c(FALSE, TRUE)  # FALSE = Scenario 1, TRUE = Scenario 2



# ==============================================================================
# Part 1: SNR sweep
# Fixed: n, M = 5, Q = 5, rho = 0.8.
# Swept: SNR in {0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}.
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_01 <- function() {
    snr_grid   <- c(0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0)
    p_p0       <- BASE$M * BASE$Q
    support_p0 <- make_support_bt_one_per_block(p = p_p0, n_blocks = BASE$M)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - SNR sweep\n")
    cat(sprintf(" n=%d p=%d M=%d Q=%d rho=%.1f SNR=swept nu=%.1f\n tFDR=%.2f beta_active=%.1f %d MC\n",
                BASE$n, p_p0, BASE$M, BASE$Q, BASE$rho, BASE$nu, BASE$tFDR, BASE$amplitude,
                BASE$num_MC))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      for (lnk in linkage_grid) {
        cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))
        snr_results <- lapply(snr_grid, function(snr_val) {
          cat(sprintf("  [SNR = %.2f] running %d MC trials ...\n", snr_val, BASE$num_MC))
          r <- .run_mc_ht_block(
            n           = BASE$n, p = p_p0, support = support_p0,
            amplitude   = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
            nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
            snr         = snr_val, tFDR = BASE$tFDR, K = BASE$K,
            num_MC      = BASE$num_MC, seed = BASE$seed,
            n_cores     = num_cores
          )
          c(list(SNR = snr_val), r)
        })
        .print_table(snr_results, "SNR", "%-6.2f")
      }
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_01()
}



# ==============================================================================
# Part 2: rho sweep
# Fixed: n, M = 5, Q = 5, SNR = 2.
# Swept: rho in {0.0, 0.1, ..., 0.9}.
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_02 <- function() {
    rho_grid <- seq(0, 0.9, by = 0.1)
    p_p1     <- BASE$M * BASE$Q
    support_p1 <- make_support_bt_one_per_block(p = p_p1, n_blocks = BASE$M)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - rho sweep\n")
    cat(sprintf(" n=%d p=%d M=%d Q=%d rho=swept SNR=%.1f nu=%.1f\n tFDR=%.2f beta_active=%.1f %d MC\n",
                BASE$n, p_p1, BASE$M, BASE$Q, BASE$snr, BASE$nu, BASE$tFDR, BASE$amplitude,
                BASE$num_MC))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      for (lnk in linkage_grid) {
        cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))
        rho_results <- lapply(rho_grid, function(rho_val) {
          cat(sprintf("  [rho = %.1f] running %d MC trials ...\n", rho_val, BASE$num_MC))
          r <- .run_mc_ht_block(
            n           = BASE$n, p = p_p1, support = support_p1,
            amplitude   = BASE$amplitude, n_blocks = BASE$M, rho = rho_val,
            nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
            snr         = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
            num_MC      = BASE$num_MC, seed = BASE$seed,
            n_cores     = num_cores
          )
          c(list(rho = rho_val), r)
        })
        .print_table(rho_results, "rho", "%-6.1f")
      }
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_02()
}



# ==============================================================================
# Part 3: Q sweep (block size / group size)
# Fixed: n, M = 5, rho = 0.8, SNR = 2.
# Swept: Q in {5, 10, 15, 20, 25, 30, 35, 40, 45, 50}.
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_03 <- function() {
    q_grid <- seq(5L, 50L, by = 5L)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - Q sweep\n")
    cat(sprintf(
      " n=%d p=swept M=%d Q=swept rho=%.1f SNR=%.1f nu=%.1f\n tFDR=%.2f beta_active=%.1f %d MC\n",
      BASE$n, BASE$M, BASE$rho, BASE$snr, BASE$nu, BASE$tFDR, BASE$amplitude, BASE$num_MC
    ))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      for (lnk in linkage_grid) {
        cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))
        q_results <- lapply(q_grid, function(q_val) {
          p_val   <- BASE$M * q_val
          support <- make_support_bt_one_per_block(p = p_val, n_blocks = BASE$M)
          cat(sprintf("  [Q = %d, p = %d] running %d MC trials ...\n",
                      q_val, p_val, BASE$num_MC))
          r <- .run_mc_ht_block(
            n           = BASE$n, p = p_val, support = support,
            amplitude   = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
            nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
            snr         = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
            num_MC      = BASE$num_MC, seed = BASE$seed,
            n_cores     = num_cores
          )
          c(list(Q = q_val, p = p_val), r)
        })
        .print_table_multi(q_results, c("Q", "p"), c(Q = "%-4d", p = "%-6d"))
      }
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_03()
}



# ==============================================================================
# Part 4: M sweep (number of groups / blocks)
# Fixed: n, Q = 5, rho = 0.8, SNR = 2.
# Swept: M in {1, 3, 5, 10, 15, 20, 30}.
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_04 <- function() {
    m_grid <- c(1L, 3L, 5L, 10L, 15L, 20L, 30L)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - M sweep\n")
    cat(sprintf(
      " n=%d p=swept M=swept Q=%d rho=%.1f SNR=%.1f nu=%.1f\n tFDR=%.2f beta_active=%.1f %d MC\n",
      BASE$n, BASE$Q, BASE$rho, BASE$snr, BASE$nu, BASE$tFDR, BASE$amplitude, BASE$num_MC
    ))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      for (lnk in linkage_grid) {
        cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))
        m_results <- lapply(m_grid, function(m_val) {
          p_val   <- m_val * BASE$Q
          support <- make_support_bt_one_per_block(p = p_val, n_blocks = m_val)
          cat(sprintf("  [M = %d, p = %d, s = %d] running %d MC trials ...\n",
                      m_val, p_val, m_val, BASE$num_MC))
          r <- .run_mc_ht_block(
            n           = BASE$n, p = p_val, support = support,
            amplitude   = BASE$amplitude, n_blocks = m_val, rho = BASE$rho,
            nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
            snr         = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
            num_MC      = BASE$num_MC, seed = BASE$seed,
            n_cores     = num_cores
          )
          c(list(M = m_val, p = p_val, s = m_val), r)
        })
        .print_table_multi(m_results, c("M", "p", "s"),
                           c(M = "%-4d", p = "%-6d", s = "%-4d"))
      }
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_04()
}



# ==============================================================================
# Part 5: Target FDR (tFDR / alpha) sweep
# Fixed: n, M = 5, Q = 5, rho = 0.8, SNR = 2.
# Swept: tFDR in seq(0.05, 0.50, by = 0.05)  (i.e., 5% to 50% in steps of 5%)
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_05 <- function() {
    tfdr_grid <- seq(0.05, 0.50, by = 0.05)
    p_p5      <- BASE$M * BASE$Q
    support_p5 <- make_support_bt_one_per_block(p = p_p5, n_blocks = BASE$M)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - tFDR sweep\n")
    cat(sprintf(
      " n=%d p=%d M=%d Q=%d rho=%.1f SNR=%.1f nu=%.1f\n tFDR=swept beta_active=%.1f %d MC\n",
      BASE$n, p_p5, BASE$M, BASE$Q, BASE$rho, BASE$snr, BASE$nu, BASE$amplitude, BASE$num_MC
    ))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      for (lnk in linkage_grid) {
        cat(sprintf("\n  --- hc_dist = \"%s\" ---\n", lnk))
        tfdr_results <- lapply(tfdr_grid, function(tfdr_val) {
          cat(sprintf("  [tFDR = %.2f] running %d MC trials ...\n", tfdr_val, BASE$num_MC))
          r <- .run_mc_ht_block(
            n           = BASE$n, p = p_p5, support = support_p5,
            amplitude   = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
            nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
            snr         = BASE$snr, tFDR = tfdr_val, K = BASE$K,
            num_MC      = BASE$num_MC, seed = BASE$seed,
            n_cores     = num_cores
          )
          c(list(tFDR = tfdr_val), r)
        })
        .print_table(tfdr_results, "tFDR", "%-6.2f")
      }
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_05()
}



# ==============================================================================
# Part 6: Linkage (Dendrogram Type) sweep
# Fixed: n, M = 5, Q = 5, rho = 0.8, SNR = 2, tFDR = 0.2
# Swept: hc_dist in {"single", "complete", "average"}
# ==============================================================================
if (TRUE) {
  trx_da_07_bt_ht_sweeps_06 <- function() {
    p_p6       <- BASE$M * BASE$Q
    support_p6 <- make_support_bt_one_per_block(p = p_p6, n_blocks = BASE$M)

    cat(strrep("=", 70), "")
    cat("\nHEAVY-TAILED X - Linkage sweep\n")
    cat(sprintf(
      " n=%d p=%d M=%d Q=%d rho=%.1f SNR=%.1f nu=%.1f\n tFDR=%.2f beta_active=%.1f %d MC\n",
      BASE$n, p_p6, BASE$M, BASE$Q, BASE$rho, BASE$snr, BASE$nu, BASE$tFDR, BASE$amplitude,
      BASE$num_MC
    ))
    cat(strrep("=", 70), "")

    for (h_noise in scenario_grid) {
      scen_str <- ifelse(h_noise, "Scenario 2: Heavy Noise", "Scenario 1: Gaussian Noise")
      cat(sprintf("\n*** %s ***\n", scen_str))

      linkage_results <- lapply(linkage_grid, function(lnk) {
        cat(sprintf("  [Linkage = %s] running %d MC trials ...\n", lnk, BASE$num_MC))
        r <- .run_mc_ht_block(
          n           = BASE$n, p = p_p6, support = support_p6,
          amplitude   = BASE$amplitude, n_blocks = BASE$M, rho = BASE$rho,
          nu          = BASE$nu, heavy_noise = h_noise, hc_dist = lnk,
          snr         = BASE$snr, tFDR = BASE$tFDR, K = BASE$K,
          num_MC      = BASE$num_MC, seed = BASE$seed,
          n_cores     = num_cores
        )
        c(list(Linkage = lnk), r)
      })
      .print_table(linkage_results, "Linkage", "%-10s")
    }
    cat("\n")
  }
  trx_da_07_bt_ht_sweeps_06()
}


cat("\nHeavy-Tailed parameter sweeps complete.\n")
cat(sprintf("\n[Log file written: %s]\n", log_file))
