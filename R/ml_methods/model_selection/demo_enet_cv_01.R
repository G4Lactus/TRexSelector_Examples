# ==============================================================================
# demo_enet_cv_01.R
# ==============================================================================
#
# Elastic-net path fitting (ElasticNet) and K-fold cross-validation
# (ElasticNetCV), both glmnet-style pathwise coordinate descent (CCD).
#
# Demo content (mirrors cpp demo_mlm_ms_02_enet_cv_ccd):
#   Part A1  ElasticNet path, low-dim  (n=300, p=100), alpha = 0 / 0.5 / 1.
#   Part A2  ElasticNet path, high-dim (n=200, p=500), alpha = 0 / 1.
#   Part B1  ElasticNetCV lambda.min / lambda.1se, low-dim, 10-fold.
#   Part B2  Ridge CV (alpha=0) + TRexGVS lambda_2_lars = lambda.1se * p / 2.
#
# Same n / p / alpha / seeds as the C++ demo; R's RNG stream differs from
# std::mt19937, so results match qualitatively, not bitwise.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("ElasticNet / ElasticNetCV Demo Suite (CCD path fit + K-fold CV)\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# Shared helpers
# ==============================================================================

# Synthetic regression data with known sparse signal: n_active coefficients
# equal to 1.0 at evenly spaced column indices, y = X beta_true + N(0, sigma^2).
# Mirrors the C++ generate_data() helper.
generate_data <- function(n, p, sigma, seed, n_active = 10L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- numeric(p)
  stride <- max(1L, p %/% n_active)
  beta_true[1L + (seq_len(n_active) - 1L) * stride] <- 1.0
  y <- as.numeric(X %*% beta_true) + sigma * rnorm(n)
  list(X = X, y = y, beta_true = beta_true)
}

# Condensed coefficient path from an ElasticNet fit: every step_size-th lambda
# with non-zero count (l0), L1-norm, and deviance ratio (glmnet's dev.ratio,
# rising as lambda decreases).
print_path_summary <- function(en, step_size = 10L) {
  lam <- en$get_lambdas()
  B   <- en$get_coef()           # p x K coefficient path
  dev <- en$get_dev_ratio()
  K   <- length(lam)

  cat(sprintf("  Path: %d lambdas fitted%s\n",
              K, if (en$converged()) "" else " [non-converged]"))
  cat(sprintf("  %-12s%-8s%-12s%-10s\n", "lambda", "nnz", "L1(beta)", "dev.ratio"))

  print_row <- function(k) {
    cat(sprintf("  %-12.5f%-8d%-12.5f%-10.5f\n",
                lam[k], sum(abs(B[, k]) > 1e-9), sum(abs(B[, k])), dev[k]))
  }
  for (k in seq(1L, K, by = step_size)) print_row(k)
  if ((K - 1L) %% step_size != 0L) print_row(K)   # always show the last point
  cat("\n")
}

# Fit an auto-grid path per alpha, print the summary, and check that (a) the
# deviance ratio is non-decreasing and (b) fit_grid() at the auto grid
# reproduces the auto path exactly.
fit_path_suite <- function(title, n, p, seed, alphas, sigma = 1.0) {
  cat(strrep("-", 70), "\n")
  cat(title, "\n")
  cat(strrep("-", 70), "\n\n")

  dat <- generate_data(n, p, sigma, seed)
  cat(sprintf("Data: n=%d, p=%d, active coefs=%d\n\n",
              n, p, sum(dat$beta_true != 0)))

  for (alpha in alphas) {
    # alpha=1 pure L1 lasso (sparse), alpha=0 pure L2 ridge (dense), in between
    # the elastic net.
    cat(sprintf("--- alpha=%g ---\n", alpha))
    en <- ElasticNet$new()$fit(dat$X, dat$y, alpha)
    print_path_summary(en, step_size = 10L)

    stopifnot(all(diff(en$get_dev_ratio()) >= -1e-12))
    en_grid  <- ElasticNet$new()$fit_grid(dat$X, dat$y, en$get_lambdas(), alpha)
    max_diff <- max(abs(en$get_coef() - en_grid$get_coef()))
    stopifnot(max_diff < 1e-12)
    cat(sprintf("  fit_grid() reproduction check (max |coef diff| = %.2e): PASS\n\n",
                max_diff))
  }
}

# Condensed summary of a fitted ElasticNetCV object:
#   - lambda.min: grid value minimizing the CV mean-squared error
#   - lambda.1se: largest lambda whose CV-MSE is within one standard error
#     of the minimum (more parsimonious / more heavily regularized choice)
# The lambda grid is stored in descending order (glmnet convention); the
# curve is printed every `step_size` grid points with the min/1se positions
# marked.
print_cv_summary <- function(cv, step_size = 10L) {
  lambdas <- cv$get_lambdas()
  mse     <- cv$get_cv_errors()
  sem     <- cv$get_cv_std()
  K       <- length(lambdas)
  i_min   <- cv$index_min()
  i_1se   <- cv$index_1se()

  cat(sprintf("  lambda.min : %.5f  (CV-MSE=%.5f +/- %.5f)\n",
              cv$cv_min(), mse[i_min], sem[i_min]))
  cat(sprintf("  lambda.1se : %.5f  (CV-MSE=%.5f +/- %.5f)\n\n",
              cv$cv_1se(), mse[i_1se], sem[i_1se]))

  cat(sprintf("  CV-MSE curve (every %d steps, lambda descending):\n",
              step_size))
  cat(sprintf("    %-14s%-12s%-12s%s\n", "lambda", "CV-MSE", "SEM", "marker"))

  print_row <- function(k, marker) {
    cat(sprintf("    %-14.5f%-12.5f%-12.5f%s\n",
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

# ==============================================================================
# Part A1: ElasticNet path — low-dimensional (n=300, p=100, n > p)
# ==============================================================================

fit_path_suite("Part A1: ElasticNet path | Low-Dimensional (n=300, p=100)",
               n = 300L, p = 100L, seed = 101L, alphas = c(0.0, 0.5, 1.0))

# ==============================================================================
# Part A2: ElasticNet path — high-dimensional (n=200, p=500, p > n)
# ==============================================================================

fit_path_suite("Part A2: ElasticNet path | High-Dimensional (n=200, p=500)",
               n = 200L, p = 500L, seed = 202L, alphas = c(0.0, 1.0))

# ==============================================================================
# Part B1: ElasticNetCV — low-dimensional CV (n=300, p=100, 10-fold)
# ==============================================================================

cat(strrep("-", 70), "\n")
cat("Part B1: ElasticNetCV | Low-Dimensional CV (n=300, p=100, 10-fold)\n")
cat(strrep("-", 70), "\n\n")

n <- 300L; p <- 100L; sigma <- 1.0; n_folds <- 10L
dat <- generate_data(n, p, sigma, seed = 303L)
cat(sprintf("Data: n=%d, p=%d, folds=%d, active coefs=%d\n\n",
            n, p, n_folds, sum(dat$beta_true != 0)))

for (alpha in c(0.0, 0.5, 1.0)) {
  # alpha mixes the penalties: alpha=1 is the pure L1 lasso (sparse fits),
  # alpha=0 pure L2 ridge (dense shrinkage), 0 < alpha < 1 the elastic net
  # in between.
  cat(sprintf("--- alpha=%g ---\n", alpha))
  cv <- ElasticNetCV$new()
  cv$fit(dat$X, dat$y, alpha = alpha, n_folds = n_folds, seed = 0L)
  print_cv_summary(cv, step_size = 10L)

  # The CV curve is evaluated at every grid point (the grid may hold fewer
  # than n_lambda points when glmnet's early termination kicks in).
  stopifnot(length(cv$get_cv_errors()) == length(cv$get_lambdas()))
  stopifnot(length(cv$get_cv_std())    == length(cv$get_lambdas()))
  stopifnot(cv$cv_1se() >= cv$cv_min())
  cat("  Ordering check (lambda.min <= lambda.1se): PASS\n\n")
}

# ==============================================================================
# Part B2: ElasticNetCV — ridge CV + lambda_2_lars conversion (n=300, p=200)
# ==============================================================================
# TRexGVS converts the glmnet-scale ridge lambda via:
#   lambda_2_lars = lambda.1se * p / 2
# mirroring R's `lm_dummy` and `cv.glmnet(alpha=0)`.
# ==============================================================================

cat(strrep("-", 70), "\n")
cat("Part B2: ElasticNetCV | Ridge CV + lambda_2_lars Conversion\n")
cat(strrep("-", 70), "\n\n")

n <- 300L; p <- 200L; sigma <- 1.0; n_folds <- 10L
dat <- generate_data(n, p, sigma, seed = 404L)
cat(sprintf("Data: n=%d, p=%d, folds=%d, active coefs=%d\n\n",
            n, p, n_folds, sum(dat$beta_true != 0)))

cv <- ElasticNetCV$new()
cv$fit(dat$X, dat$y, alpha = 0.0, n_folds = n_folds, seed = 0L)
print_cv_summary(cv, step_size = 10L)

stopifnot(cv$cv_1se() >= cv$cv_min())
lambda2_min <- cv$cv_min() * p / 2
lambda2_1se <- cv$cv_1se() * p / 2

cat(sprintf("  lambda_2_lars (min) = cv_min  * p/2 = %.6f * %d/2 = %.6f\n",
            cv$cv_min(), p, lambda2_min))
cat(sprintf("  lambda_2_lars (1se) = cv_1se  * p/2 = %.6f * %d/2 = %.6f\n\n",
            cv$cv_1se(), p, lambda2_1se))
cat("  (This is the value TRexGVS uses internally via CV_1SE_CCD.)\n\n")

cat(strrep("=", 70), "\n")
cat("demo_enet_cv_01.R complete.\n")
