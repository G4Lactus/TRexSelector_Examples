# ==============================================================================
# demo_trex_spca_sim_parallel.R
# ==============================================================================
#
# This script evaluates the T-Rex Sparse PCA method against ordinary PCA and an
# oracle baseline.
#
# 1. Data Generating Process (DGP)
#    - Model: Sparse M-Factor Model (X = Z * V^T + E)
#    - Dimensions: n = 50 samples, p = 100 features.
#    - Factors (Z): M = 3 independent latent factors with decreasing
#      standard deviations (sigma = 5, 3, 1).
#    - Sparse Loadings (V): Each factor has p1 = 5 true active variables.
#      To enforce overlapping support, these variables are strictly sampled
#      from the same pool of the first 30 variables.
#    - Noise (E): Scaled according to Signal-to-Noise Ratio (SNR) targets
#      ranging from -10 dB to +10 dB.
#
# 2. Evaluation Strategy
#    - TPR & FDR (Support Recovery): Because the true factors share overlapping
#      variables, Ordinary PCA (enforces strict orthogonality) forces
#      the 2nd and 3rd estimated principal components to become mixed linear
#      combinations of the true factors.
#      To avoid penalizing algorithms for finding these "leaked" variables,
#      TPR and FDR are strictly evaluated on the 1st PC (PC1)
#      where the dominant factor (sigma = 5) isolates cleanly.
#    - Percentage of Explained Variance (PEV): Evaluated cumulatively across
#      all M components. Because sparse principal components lose strict
#      orthogonality, PEV is adjusted using the QR decomposition of the
#      estimated score matrix (Z_est) as defined in the paper.
#
# 3. Paper Visualizations
#    While this script sweeps across SNR values to output console metrics,
#    the original paper plots these relationships across multiple axes:
#      * Fig 1/2: TPR & FDR vs. SNR (dB)
#      * Fig 3a:  Cumulative PEV vs. SNR (dB)
#      * Fig 3b:  Cumulative PEV vs. Number of true active loadings (p1)
#      * Fig 3c:  Cumulative PEV vs. Target FDR (alpha)
#
#
# NOTE ON T-REX+GVS with ELASTIC NET PARAMETERIZATION:
# --------------------------------------------------------
# To recover the highly correlated overlapping variables without
# reverting to Lasso-style selection (picking 1 out of 5), the internal
# selector is the TRex+GVS selector operating an Elastic Net base selector.
# We rely on the internal cross-validation (lambda_2_lars = NULL) to
# autonomously capture the grouping effect.
# ==============================================================================

library(TRexSelector)
library(elasticnet)
library(parallel)

# ------------------------------------------------------------------------------
# 1. Data Generating Process (Sparse M-Factor Model)
# ------------------------------------------------------------------------------
dgp_sparse_factor_model <- function(n = 50, p = 100, p1 = 5, M = 3,
                                    target_snr_db = 0,
                                    overlap_pool = 30) {

  Z <- matrix(0, nrow = n, ncol = M)
  factor_stds <- c(5.0, 3.0, 1.0)
  for (m in 1:M) {
    std_dev <- ifelse(m <= length(factor_stds), factor_stds[m], 1.0)
    Z[, m] <- rnorm(n, mean = 0, sd = std_dev)
  }

  V <- matrix(0, nrow = p, ncol = M)
  pool_indices <- 1:min(overlap_pool, p)
  for (m in 1:M) {
    active_idx <- sample(pool_indices, p1, replace = FALSE)
    V[active_idx, m] <- 0.9
  }

  S <- Z %*% t(V)
  signal_var <- var(as.vector(S))

  E_raw <- matrix(rnorm(n * p), nrow = n, ncol = p)
  raw_noise_var <- var(as.vector(E_raw))

  target_noise_var <- signal_var / (10^(target_snr_db / 10))
  E <- E_raw * sqrt(target_noise_var / raw_noise_var)

  X <- S + E
  X <- scale(X, center = TRUE, scale = FALSE)
  list(X = X, V = V, Z = Z)
}

# ------------------------------------------------------------------------------
# 2. Utils Evaluator
# ------------------------------------------------------------------------------
evaluate_spca <- function(X, V_est, Z_est, V_true) {
  n <- nrow(X)
  p <- ncol(X)
  M <- ncol(V_est)

  true_supp_pc1 <- which(V_true[, 1] != 0)
  est_supp_pc1  <- which(abs(V_est[, 1]) > 1e-12)

  tp <- length(intersect(true_supp_pc1, est_supp_pc1))
  tpr <- ifelse(length(true_supp_pc1) == 0, 0.0,
                tp / max(1, length(true_supp_pc1)))
  fdr <- ifelse(length(est_supp_pc1) == 0, 0.0,
                (length(est_supp_pc1) - tp) / max(1, length(est_supp_pc1)))

  qr_Z <- qr(Z_est)
  R <- qr.R(qr_Z)

  total_var <- sum(X^2) / (n - 1)
  cum_pev <- numeric(M)
  cum_ev <- 0

  for (m in 1:M) {
    r_mm <- R[m, m]
    cum_ev <- cum_ev + (r_mm^2) / (n - 1)
    cum_pev[m] <- cum_ev / total_var
  }

  list(TPR = tpr, FDR = fdr, PEV = cum_pev)
}

# ------------------------------------------------------------------------------
# 3. Main Monte Carlo Simulation (Parallelized)
# ------------------------------------------------------------------------------
n <- 50
p <- 100
p1 <- 5
M <- 3
target_fdr <- 0.10
num_MC <- 200
snr_db_seq <- c(-10, -7, -5, -3, 0, 3, 5, 7, 10)

num_cores <- 8
is_unix <- .Platform$OS.type == "unix"

cat("========================================================\n")
cat(sprintf("  T-Rex SPCA MC Simulation (%d cores)\n", num_cores))
cat("========================================================\n\n")

# Initialize Windows Cluster once to avoid overhead
if (!is_unix) {
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, {
    library(TRexSelector)
    library(elasticnet)
  })
  clusterExport(cl, varlist = c("n", "p", "p1", "M", "target_fdr",
                                "dgp_sparse_factor_model", "evaluate_spca"),
                envir = environment())
}

for (snr in snr_db_seq) {
  cat(sprintf("Running SNR = %3d dB...\n", snr))

  if (!is_unix) {
    clusterExport(cl, "snr", envir = environment())
  }

  one_trial <- function(mc) {
    set.seed(42 + mc * 1000 + snr)
    data <- dgp_sparse_factor_model(n, p, p1, M, snr)
    X <- data$X
    V_true <- data$V

    # Expanded to 24 elements to accommodate 8 distinct methods
    res <- numeric(24)
    names(res) <- c(
      "ord_TPR", "ord_FDR", "ord_PEV",
      "oracle_TPR", "oracle_FDR", "oracle_PEV",
      "en_spca_TPR",  "en_spca_FDR",  "en_spca_PEV",
      "oracle_spca_TPR", "oracle_spca_FDR", "oracle_spca_PEV",
      "en_act_TPR", "en_act_FDR", "en_act_PEV",
      "en_thr_TPR", "en_thr_FDR", "en_thr_PEV"
    )

    # ---------------------------------------------------------
    # Method 1: Ordinary PCA
    # ---------------------------------------------------------
    svd_X <- svd(X)
    V_ord <- svd_X$v[, 1:M]
    Z_ord <- X %*% V_ord

    ev <- evaluate_spca(X, V_ord, Z_ord, V_true)
    res["ord_TPR"] <- ev$TPR
    res["ord_FDR"] <- ev$FDR
    res["ord_PEV"] <- ev$PEV[M]


    # ---------------------------------------------------------
    # Method 2: Oracle Thresholded PCA
    # ---------------------------------------------------------
    V_oracle <- matrix(0, nrow = p, ncol = M)
    for (m in 1:M) {
      v_m <- V_ord[, m]
      thresh <- sort(abs(v_m), decreasing = TRUE)[p1]
      v_m[abs(v_m) < thresh] <- 0
      V_oracle[, m] <- v_m / sqrt(sum(v_m^2))
    }
    ev <- evaluate_spca(X, V_oracle, X %*% V_oracle, V_true)
    res["oracle_TPR"] <- ev$TPR
    res["oracle_FDR"] <- ev$FDR
    res["oracle_PEV"] <- ev$PEV[M]


    # ---------------------------------------------------------
    # Method 3: Zou et al. (2006) Elastic Net SPCA
    # ---------------------------------------------------------
    en_res <- elasticnet::spca(X, K = M, para = rep(p1, M),
                               type = "predictor",
                               sparse = "varnum",
                               trace = FALSE)
    V_en_spca <- en_res$loadings
    Z_en_spca <- X %*% V_en_spca

    ev_en_spca <- evaluate_spca(X, V_en_spca, Z_en_spca, V_true)
    res["en_spca_TPR"] <- ev_en_spca$TPR
    res["en_spca_FDR"] <- ev_en_spca$FDR
    res["en_spca_PEV"] <- ev_en_spca$PEV[M]


    # ---------------------------------------------------------
    # Method 4: Oracle SPCA
    # ---------------------------------------------------------
    V_oracle_spca <- matrix(0, p, M)
    for (m in 1:M) {
      y_m <- Z_ord[, m]

      # 1. Force heavy Ridge penalty to activate grouping
      en_fit <- elasticnet::enet(x = X, y = y_m,
                                 lambda = 1000,
                                 normalize = FALSE, trace = FALSE)

      # 2. First step where p1 variables are active
      active_counts <- apply(en_fit$beta, 1, function(b) sum(abs(b) > 1e-12))
      step_idx <- which(active_counts >= p1)[1]

      # Fallback to the end of the path if they never reach p1
      if (is.na(step_idx)) step_idx <- length(active_counts)

      # 3. Extract the physical coefficients
      beta_en <- en_fit$beta[step_idx, ]

      # 4. Strict Cardinality Enforcement: Prune down to exactly p1 variables
      thresh <- sort(abs(beta_en), decreasing = TRUE)[p1]
      beta_en[abs(beta_en) < thresh] <- 0

      if (sum(beta_en^2) > 0) {
        V_oracle_spca[, m] <- beta_en / sqrt(sum(beta_en^2))
      }
    }

    ev_oracle_spca <- evaluate_spca(X, V_oracle_spca, X %*% V_oracle_spca, V_true)
    res["oracle_spca_TPR"] <- ev_oracle_spca$TPR
    res["oracle_spca_FDR"] <- ev_oracle_spca$FDR
    res["oracle_spca_PEV"] <- ev_oracle_spca$PEV[M]


    # ---------------------------------------------------------
    # Method 5 & 6: T-Rex SPCA (Elastic Net GVS)
    # ---------------------------------------------------------
    V_en_active <- matrix(0, p, M)
    V_en_thresh <- matrix(0, p, M)

    for (m in 1:M) {
      y_m <- Z_ord[, m]

      trex_res_en <- TRexSelector::trex(X = X, y = y_m, tFDR = target_fdr,
                                        method = "trex+GVS", GVS_type = "EN",
                                        lambda_2_lars = NULL, verbose = FALSE)

      active_set <- which(trex_res_en$selected_var > .Machine$double.eps)
      k <- length(active_set)

      if (k > 0) {
        v_m <- V_ord[, m]
        thresh <- sort(abs(v_m), decreasing = TRUE)[k]
        v_m[abs(v_m) < thresh] <- 0
        V_en_thresh[, m] <- v_m / sqrt(sum(v_m^2))

        X_A <- X[, active_set, drop = FALSE]
        lambda2 <- 1e-6
        I_k <- diag(1, nrow = k, ncol = k)
        beta_ridge <- solve(crossprod(X_A) + lambda2 * I_k) %*% crossprod(X_A, y_m)
        V_en_active[active_set, m] <- beta_ridge / sqrt(sum(beta_ridge^2))
      }
    }

    ev_en_a <- evaluate_spca(X, V_en_active, X %*% V_en_active, V_true)
    res["en_act_TPR"] <- ev_en_a$TPR
    res["en_act_FDR"] <- ev_en_a$FDR
    res["en_act_PEV"] <- ev_en_a$PEV[M]

    ev_en_t <- evaluate_spca(X, V_en_thresh, X %*% V_en_thresh, V_true)
    res["en_thr_TPR"] <- ev_en_t$TPR
    res["en_thr_FDR"] <- ev_en_t$FDR
    res["en_thr_PEV"] <- ev_en_t$PEV[M]


    return(res)
  }

  # Execute parallel MC trials
  if (is_unix) {
    res_list <- mclapply(seq_len(num_MC), one_trial, mc.cores = num_cores, mc.preschedule = FALSE)
  } else {
    res_list <- parLapply(cl, seq_len(num_MC), one_trial)
  }


  # Aggregate results
  # -----------------------------------
  res_mat <- do.call(rbind, res_list)
  avg <- colMeans(res_mat)

  # Print results
  # -----------------------------------
  cat(sprintf("  -> Ord PCA        | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["ord_TPR"] * 100, avg["ord_FDR"] * 100, avg["ord_PEV"] * 100))
  cat(sprintf("  -> Oracle PCA     | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["oracle_TPR"] * 100, avg["oracle_FDR"] * 100, avg["oracle_PEV"] * 100))
  cat(sprintf("  -> Zou EN SPCA    | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["en_spca_TPR"] * 100, avg["en_spca_FDR"] * 100, avg["en_spca_PEV"] * 100))
  cat(sprintf("  -> Oracle SPCA    | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["oracle_spca_TPR"] * 100, avg["oracle_spca_FDR"] * 100, avg["oracle_spca_PEV"] * 100))
  cat(sprintf("  -> T-Rex EN Act   | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["en_act_TPR"] * 100, avg["en_act_FDR"] * 100, avg["en_act_PEV"] * 100))
  cat(sprintf("  -> T-Rex EN Thr   | TPR: %5.1f%% | FDR: %5.1f%% | PEV(M=3): %5.1f%%\n",
              avg["en_thr_TPR"] * 100, avg["en_thr_FDR"] * 100, avg["en_thr_PEV"] * 100))
}
