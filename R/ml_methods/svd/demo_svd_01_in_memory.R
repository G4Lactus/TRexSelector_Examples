# ==============================================================================
# demo_svd_01_in_memory.R
# ==============================================================================
#
# Demo illustrating SVDSolver (truncated rank-M SVD) on in-memory data.
#
# Demo content:
#   - Part 1: Direct path (n >= p) — truncated and full-rank reconstruction
#   - Part 2: Gram / kernel path (p > 2*n) — wide-matrix decomposition
#   - Part 3: Orthogonality check — ||U^T U - I||_F and ||V^T V - I||_F
#
# Mirrors: cpp/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01.cpp
#
# Note: same n / p / M values and seeds as the C++ demo; R's Mersenne-Twister
# stream differs from std::mt19937, so results match qualitatively, not
# bitwise.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("SVDSolver Demos (in-memory)\n")
cat(strrep("=", 70), "\n\n")

# The truncated SVD keeps the M leading singular triplets of X (n x p):
#   X_M = U_M * diag(S_M) * V_M^T,   U_M: n x M,  S_M: M,  V_M: p x M
# X_M is the best rank-M approximation of X in the Frobenius norm
# (Eckart-Young). $compute(X, M) returns a list with elements U, S, V.

solver <- SVDSolver$new()

# ==============================================================================
# Part 1: Direct / Thin Path (n >= p)
# ==============================================================================
# For tall matrices the solver decomposes X directly.
# ==============================================================================

cat(strrep("-", 70), "\n")
cat("Part 1: Direct Path (n >= p)\n")
cat(strrep("-", 70), "\n\n")

set.seed(42)
n1 <- 300L
p1 <- 80L
M1 <- 10L

cat(sprintf("Matrix size : %d x %d\n", n1, p1))
cat(sprintf("M           : %d\n\n", M1))

X1 <- matrix(rnorm(n1 * p1), nrow = n1, ncol = p1)

res1 <- solver$compute(X1, M1)

cat(sprintf("U dimensions : %d x %d\n", nrow(res1$U), ncol(res1$U)))
cat(sprintf("S length     : %d\n", length(res1$S)))
cat(sprintf("V dimensions : %d x %d\n\n", nrow(res1$V), ncol(res1$V)))

# Relative reconstruction error of the rank-M approximation:
# ||X_M - X||_F / ||X||_F  (large here: a Gaussian matrix has a flat
# spectrum, so 10 of 80 components capture little energy).
X1_approx <- res1$U %*% diag(res1$S) %*% t(res1$V)
rel_error1 <- norm(X1_approx - X1, "F") / norm(X1, "F")
cat("Relative reconstruction error ||X_M - X||_F / ||X||_F :",
    format(rel_error1, scientific = TRUE), "\n")

cat(sprintf("\nTop %d singular values:\n", M1))
for (k in seq_len(M1)) {
  cat(sprintf("  s%d = %.4f\n", k, res1$S[k]))
}

# Didactic checks: singular values are sorted decreasingly, and at full rank
# (M = p) the truncated SVD reproduces X exactly up to floating-point error.
stopifnot(all(diff(res1$S) <= 1e-12))

res1_full <- solver$compute(X1, p1)
X1_full   <- res1_full$U %*% diag(res1_full$S) %*% t(res1_full$V)
rel_error_full <- norm(X1_full - X1, "F") / norm(X1, "F")
cat("\nFull-rank (M = p) relative reconstruction error :",
    format(rel_error_full, scientific = TRUE), "\n")
stopifnot(rel_error_full < 1e-10)
cat("Check passed: full-rank reconstruction is exact.\n")

# ==============================================================================
# Part 2: Gram / Kernel Path (p > 2*n)
# ==============================================================================
# For wide matrices the solver works on the n x n Gram matrix X X^T instead
# of X itself, which is much cheaper when p >> n.
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 2: Gram Path (p > 2*n)\n")
cat(strrep("-", 70), "\n\n")

set.seed(7)
n2 <- 100L
p2 <- 500L
M2 <- 8L

cat(sprintf("Matrix size : %d x %d  (p > 2*n -> Gram path)\n", n2, p2))
cat(sprintf("M           : %d\n\n", M2))

X2 <- matrix(rnorm(n2 * p2), nrow = n2, ncol = p2)

res2 <- solver$compute(X2, M2)

cat(sprintf("U dimensions : %d x %d\n", nrow(res2$U), ncol(res2$U)))
cat(sprintf("S length     : %d\n", length(res2$S)))
cat(sprintf("V dimensions : %d x %d\n\n", nrow(res2$V), ncol(res2$V)))

X2_approx <- res2$U %*% diag(res2$S) %*% t(res2$V)
rel_error2 <- norm(X2_approx - X2, "F") / norm(X2, "F")
cat("Relative reconstruction error ||X_M - X||_F / ||X||_F :",
    format(rel_error2, scientific = TRUE), "\n")

# ==============================================================================
# Part 3: Orthogonality Check
# ==============================================================================
# The singular vector blocks must be orthonormal: U^T U = I and V^T V = I.
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 3: Orthogonality Check\n")
cat(strrep("-", 70), "\n\n")

set.seed(99)
n3 <- 200L
p3 <- 150L
M3 <- 20L

cat(sprintf("Matrix size : %d x %d\n", n3, p3))
cat(sprintf("M           : %d\n\n", M3))

X3 <- matrix(rnorm(n3 * p3), nrow = n3, ncol = p3)

res3 <- solver$compute(X3, M3)

I_M       <- diag(M3)
UtU_error <- norm(crossprod(res3$U) - I_M, "F")
VtV_error <- norm(crossprod(res3$V) - I_M, "F")

cat("||U^T U - I||_F :", format(UtU_error, scientific = TRUE), "\n")
cat("||V^T V - I||_F :", format(VtV_error, scientific = TRUE), "\n")

stopifnot(UtU_error < 1e-10, VtV_error < 1e-10)
cat("\nCheck passed: U and V are orthonormal.\n")

cat("\n", strrep("=", 70), "\n")
cat("demo_svd_01_in_memory.R complete.\n")
