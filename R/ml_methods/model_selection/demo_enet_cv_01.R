# ==============================================================================
# demo_enet_cv_01.R
# ==============================================================================
#
# Demo illustrating elastic-net K-fold cross-validation via ElasticNetCV
# (glmnet-style pathwise cyclic coordinate descent, CCD).
#
# Demo content:
#   - Part B1: CV lambda selection (lambda.min / lambda.1se) for
#     alpha = {0.0 (ridge), 0.5 (EN), 1.0 (lasso)},
#     low-dimensional (n=300, p=100, 10-fold).
#   - Part B2: CV for pure ridge (alpha=0) as used by TRexGVS; shows the
#     lambda_2_lars = lambda.1se * p / 2 conversion.
#
# Mirrors: cpp/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/
#          demo_mlm_ms_02_enet_cv_ccd.cpp
#
# Notes on the port:
#   - Same n / p / alpha / fold values and data seeds as the C++ demo; R's
#     Mersenne-Twister stream differs from std::mt19937, so results match
#     qualitatively, not bitwise.
#   - The C++ demo additionally has a Part A (elastic-net PATH fitting via
#     `enet_gaussian` for low- and high-dimensional designs). TRexSelectorNeo's
#     R bindings expose only the CV class (ElasticNetCV) and no standalone
#     path class, so Part A is omitted here; see the Python port
#     (Python/ml_methods/model_selection/demo_enet_cv_01.py) for the path
#     scenarios.
#   - ElasticNetCV builds the full-data lambda grid once and reuses it
#     across folds (glmnet semantics). The grid is stored in glmnet's
#     DESCENDING order (largest penalty first) — unlike RidgeCV's ascending
#     grid in demo_ridge_cv_01.R — and may hold fewer than n_lambda points
#     when glmnet's early-termination rule stops the path early.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("ElasticNetCV Demo Suite (K-fold CV via coordinate descent)\n")
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
