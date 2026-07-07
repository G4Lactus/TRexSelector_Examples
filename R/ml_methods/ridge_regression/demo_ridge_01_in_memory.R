# ==============================================================================
# demo_ridge_01_in_memory.R
# ==============================================================================
#
# Demo illustrating RidgeSolver on in-memory data.
#
# Demo content:
#   - Part 1: Primal path (n >= p) — recovering a dense linear model
#   - Part 2: Dual path (n < p) — high-dimensional, sparse ground truth
#   - Part 3: Lambda sweep — coefficient shrinkage as the penalty grows
#
# Mirrors: cpp/ml_methods/ridge_regression/demo_mlm_ridge_01/demo_mlm_ridge_01.cpp
#
# Note: same n / p / lambda values and seeds as the C++ demo; R's
# Mersenne-Twister stream differs from std::mt19937, so results match
# qualitatively, not bitwise.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("RidgeSolver Demos (in-memory)\n")
cat(strrep("=", 70), "\n\n")

# Ridge regression solves the penalized least-squares problem
#   beta_hat = argmin ||y - X beta||_2^2 + lambda ||beta||_2^2 .
# The solver picks its computational path automatically:
#   - primal (p x p system) when n >= p,
#   - dual   (n x n system) when n <  p — same solution, cheaper factorization.

solver <- RidgeSolver$new()

l2_norm <- function(v) sqrt(sum(v^2))

# ==============================================================================
# Part 1: Basic Solve (n >= p, primal path)
# ==============================================================================
# Dense ground truth (all coefficients = 1), low noise: with mild
# regularization the ridge estimate should sit close to beta_true.
# ==============================================================================

cat(strrep("-", 70), "\n")
cat("Part 1: Primal Path (n >= p)\n")
cat(strrep("-", 70), "\n\n")

set.seed(42)
n1      <- 500L
p1      <- 50L
lambda1 <- 1.0

cat(sprintf("Matrix size : %d x %d\n", n1, p1))
cat(sprintf("Lambda      : %g\n\n", lambda1))

X1         <- matrix(rnorm(n1 * p1), nrow = n1, ncol = p1)
beta_true1 <- rep(1.0, p1)
y1         <- as.numeric(X1 %*% beta_true1) + 0.1 * rnorm(n1)

beta_hat1 <- solver$solve(X1, y1, lambda1)

coef_error1 <- l2_norm(beta_hat1 - beta_true1)
residual1   <- l2_norm(y1 - X1 %*% beta_hat1)
cat("||beta_hat - beta_true||_2 :", format(coef_error1, scientific = TRUE), "\n")
cat("Residual ||y - X*beta||_2  :", format(residual1, scientific = TRUE), "\n")

# Didactic check: with sd(noise) = 0.1 the estimate must be close to truth.
stopifnot(coef_error1 < 0.5)
cat("\nCheck passed: ridge estimate close to the true coefficients.\n")

# ==============================================================================
# Part 2: Dual Path (n < p)
# ==============================================================================
# High-dimensional regime: only the first 10 of 300 coefficients are nonzero.
# Ridge does not select variables — it shrinks everything — so the estimate
# is biased, but the dual formulation keeps the problem well posed for n < p.
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 2: Dual Path (n < p)\n")
cat(strrep("-", 70), "\n\n")

set.seed(7)
n2      <- 80L
p2      <- 300L
lambda2 <- 10.0

cat(sprintf("Matrix size : %d x %d\n", n2, p2))
cat(sprintf("Lambda      : %g\n\n", lambda2))

X2 <- matrix(rnorm(n2 * p2), nrow = n2, ncol = p2)

# Sparse ground truth: only the first 10 coefficients are nonzero.
beta_true2 <- c(rep(1.0, 10L), rep(0.0, p2 - 10L))
y2         <- as.numeric(X2 %*% beta_true2) + 0.1 * rnorm(n2)

beta_hat2 <- solver$solve(X2, y2, lambda2)

coef_error2 <- l2_norm(beta_hat2 - beta_true2)
residual2   <- l2_norm(y2 - X2 %*% beta_hat2)
cat("||beta_hat - beta_true||_2 :", format(coef_error2, scientific = TRUE), "\n")
cat("Residual ||y - X*beta||_2  :", format(residual2, scientific = TRUE), "\n")

# ==============================================================================
# Part 3: Lambda Sweep
# ==============================================================================
# As lambda grows, ||beta_hat||_2 shrinks monotonically toward 0 while the
# training residual grows — the classic bias-variance trade-off knob.
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 3: Lambda Sweep\n")
cat(strrep("-", 70), "\n\n")

set.seed(99)
n3 <- 200L
p3 <- 40L

cat(sprintf("Matrix size : %d x %d\n\n", n3, p3))

X3         <- matrix(rnorm(n3 * p3), nrow = n3, ncol = p3)
beta_true3 <- rep(1.0, p3)
y3         <- as.numeric(X3 %*% beta_true3) + 0.5 * rnorm(n3)

lambdas <- c(1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0)

cat(sprintf("%12s  %20s  %16s  %16s\n",
            "lambda", "||beta_hat - true||_2", "||beta_hat||_2", "||residual||_2"))
cat(strrep("-", 70), "\n")

beta_norms <- numeric(length(lambdas))
residuals  <- numeric(length(lambdas))
for (i in seq_along(lambdas)) {
  beta_hat      <- solver$solve(X3, y3, lambdas[i])
  beta_norms[i] <- l2_norm(beta_hat)
  residuals[i]  <- l2_norm(y3 - X3 %*% beta_hat)
  cat(sprintf("%12.4e  %20.6e  %16.6e  %16.6e\n",
              lambdas[i], l2_norm(beta_hat - beta_true3),
              beta_norms[i], residuals[i]))
}

# Didactic checks: for ridge, ||beta_hat(lambda)||_2 is non-increasing and
# the training residual is non-decreasing in lambda.
stopifnot(all(diff(beta_norms) <= 1e-10))
stopifnot(all(diff(residuals) >= -1e-10))
cat("\nCheck passed: coefficients shrink and residual grows with lambda.\n")

cat("\n", strrep("=", 70), "\n")
cat("demo_ridge_01_in_memory.R complete.\n")
