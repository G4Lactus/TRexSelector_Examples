# ==============================================================================
# demo_ridge_cv_01.R
# ==============================================================================
#
# Demo illustrating ridge-regression K-fold cross-validation via RidgeCV
# (SVD-based solver with per-fold centering and column L2-normalization).
#
# Demo content:
#   - Scenario 1: Low-dimensional  (n=1000, p=500,  n > p,  10-fold)
#   - Scenario 2: High-dimensional (n=300,  p=500,  p > n,  10-fold)
#   - Scenario 3: Very high-dim    (n=300,  p=1000, p >> n,  5-fold)
#   - Each scenario reports the lambda grid endpoints, the CV-MSE curve,
#     and the lambda.min / lambda.1se selections.
#
# Mirrors: cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/
#          demo_mlm_ms_01_ridge_cv_svd.cpp
#
# Notes on the port:
#   - Same n / p / fold values and data seeds as the C++ demo; R's
#     Mersenne-Twister stream differs from std::mt19937, so results match
#     qualitatively, not bitwise.
#   - The C++ demo has a 4th, memory-mapped scenario (the design matrix is
#     read from disk via a zero-copy Eigen map). RidgeCV's R binding only
#     accepts in-memory matrices, so that scenario is omitted here.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("RidgeCV Demo Suite (K-fold CV for the ridge penalty)\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# Shared helpers
# ==============================================================================

# Synthetic regression data with known sparse signal: n_active coefficients
# equal to 1.0 at evenly spaced column indices, y = X beta_true + N(0, sigma^2).
# Mirrors the C++ generate_regression_data() helper.
generate_regression_data <- function(n, p, sigma, seed, n_active = 10L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- numeric(p)
  stride <- max(1L, p %/% n_active)
  beta_true[1L + (seq_len(n_active) - 1L) * stride] <- 1.0
  y <- as.numeric(X %*% beta_true) + sigma * rnorm(n)
  list(X = X, y = y, beta_true = beta_true)
}

# Condensed summary of a fitted RidgeCV object:
#   - lambda.min: grid value minimizing the CV mean-squared error
#   - lambda.1se: largest lambda whose CV-MSE is within one standard error
#     of the minimum (more parsimonious / more heavily regularized choice)
# The lambda grid is stored in ascending order; the curve is printed every
# `step_size` grid points with the min/1se positions marked.
print_cv_summary <- function(cv, step_size = 10L) {
  lambdas <- cv$get_lambdas()
  mse     <- cv$get_cv_errors()
  sem     <- cv$get_cv_std()
  K       <- length(lambdas)
  i_min   <- cv$index_min()
  i_1se   <- cv$index_1se()

  cat(sprintf("  lambda.min : %.4f  (CV-MSE=%.4f +/- %.4f)\n",
              cv$cv_min(), mse[i_min], sem[i_min]))
  cat(sprintf("  lambda.1se : %.4f  (CV-MSE=%.4f +/- %.4f)\n\n",
              cv$cv_1se(), mse[i_1se], sem[i_1se]))

  cat(sprintf("  CV-MSE curve (every %d steps, lambda ascending):\n", step_size))
  cat(sprintf("    %-14s%-12s%-12s%s\n", "lambda", "CV-MSE", "SEM", "marker"))

  print_row <- function(k, marker) {
    cat(sprintf("    %-14.4f%-12.4f%-12.4f%s\n",
                lambdas[k], mse[k], sem[k], marker))
  }
  shown <- seq(1L, K, by = step_size)
  for (k in shown) {
    marker <- if (k == i_min) "  <-- min" else if (k == i_1se) "  <-- 1se" else ""
    print_row(k, marker)
  }
  if (!(i_min %in% shown)) print_row(i_min, "  <-- min")
  if (!(i_1se %in% shown) && i_1se != i_min) print_row(i_1se, "  <-- 1se")
  cat("\n")
}

# Runs one CV scenario end to end and applies the ordering check
# lambda.1se >= lambda.min (the 1se rule always regularizes at least as
# strongly as the pure minimizer).
run_cv_scenario <- function(title, n, p, num_folds, seed, sigma = 1.0) {
  cat(strrep("-", 70), "\n")
  cat(title, "\n")
  cat(strrep("-", 70), "\n\n")

  dat <- generate_regression_data(n, p, sigma, seed)
  cat(sprintf("Data: n=%d, p=%d, folds=%d, active coefs=%d\n\n",
              n, p, num_folds, sum(dat$beta_true != 0)))

  # Defaults mirror the C++ fit(): 100-lambda grid anchored at the largest
  # useful penalty, spanning a 1000x range; `seed` fixes fold assignment.
  cv <- RidgeCV$new()
  cv$fit(dat$X, dat$y, num_folds = num_folds,
         n_lambda = 100L, lambda_ratio = 1000, seed = seed)

  print_cv_summary(cv, step_size = 10L)

  stopifnot(cv$cv_1se() >= cv$cv_min())
  cat("  Ordering check (lambda.min <= lambda.1se): PASS\n\n")
  invisible(cv)
}

# ==============================================================================
# Scenario 1: Low-dimensional (n=1000, p=500, n > p, 10-fold)
# ==============================================================================

run_cv_scenario("Scenario 1: Low-Dimensional (n=1000, p=500, 10-fold)",
                n = 1000L, p = 500L, num_folds = 10L, seed = 111L)

# ==============================================================================
# Scenario 2: High-dimensional (n=300, p=500, p > n, 10-fold)
# ==============================================================================

run_cv_scenario("Scenario 2: High-Dimensional (n=300, p=500, 10-fold)",
                n = 300L, p = 500L, num_folds = 10L, seed = 222L)

# ==============================================================================
# Scenario 3: Very high-dimensional (n=300, p=1000, p >> n, 5-fold)
# ==============================================================================
# 5 folds keep each fold's training set large enough relative to n.
# ==============================================================================

run_cv_scenario("Scenario 3: Very High-Dimensional (n=300, p=1000, 5-fold)",
                n = 300L, p = 1000L, num_folds = 5L, seed = 333L)

cat(strrep("=", 70), "\n")
cat("demo_ridge_cv_01.R complete.\n")
