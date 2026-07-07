# ==============================================================================
# demo_pca_01_in_memory.R
# ==============================================================================
#
# Demo illustrating Principal Component Analysis (PCA) on in-memory data.
#
# Demo content:
#   - Part 1: Basic PCA fit — score/loading dimensions and explained variance
#   - Part 2: Restore round-trip accuracy: Data -> PCA -> Restore -> Data
#   - Part 3: Transforming new (test) data into the fitted PC space
#
# Mirrors: cpp/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01.cpp
#
# Notes on the port:
#   - Same n / p / M values and seeds as the C++ demo; R's Mersenne-Twister
#     stream differs from std::mt19937, so results match qualitatively,
#     not bitwise.
#   - In-place (zero-copy) design: PCA$new(X, ...) maps the memory of the R
#     matrix X directly, the same Eigen::Map interface that serves both
#     in-memory and memory-mapped data. This is deliberate — the library
#     targets data volumes that can exceed RAM, so it never takes hidden
#     copies of X. Consequently $fit() standardizes X in place and $restore()
#     undoes the preprocessing exactly. A plain alias (`X_orig <- X`) shares
#     the same buffer; where a pristine copy is needed, take an explicit deep
#     copy with `X_orig <- X + 0`.
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("PCA Demos (in-memory)\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# Part 1: Basic PCA Fit
# ==============================================================================
# PCA projects X (n x p) onto the M orthogonal directions of maximal variance.
# After fit():
#   - scores   Z (n x M): coordinates of each sample in the PC basis
#   - loadings V (p x M): the principal directions (columns are orthonormal)
#   - explained_variance (M): variance captured by each retained component
# ==============================================================================

cat(strrep("-", 70), "\n")
cat("Part 1: Basic PCA Fit\n")
cat(strrep("-", 70), "\n\n")

set.seed(42)
n1 <- 500L
p1 <- 200L
M1 <- 10L

cat(sprintf("Matrix size  : %d x %d\n", n1, p1))
cat(sprintf("Components M : %d\n\n", M1))

X1 <- matrix(rnorm(n1 * p1), nrow = n1, ncol = p1)

# center = TRUE subtracts column means; normalize = TRUE scales columns to
# unit L2 norm. Both happen in place on X1 and are undone by restore().
pca1 <- PCA$new(X1, M = M1, center = TRUE, normalize = TRUE)
pca1$fit()

Z1 <- pca1$get_scores()
V1 <- pca1$get_loadings()
ev <- pca1$get_explained_variance()

cat(sprintf("Z dimensions : %d x %d\n", nrow(Z1), ncol(Z1)))
cat(sprintf("V dimensions : %d x %d\n", nrow(V1), ncol(V1)))
cat(sprintf("explained_variance length: %d\n\n", length(ev)))

# Explained variance as a percentage of the variance retained by the M
# components (the C++ demo reports the same "of retained" normalization).
total_var <- sum(ev)
cat("Explained variance per component (% of retained):\n")
for (k in seq_len(M1)) {
  cat(sprintf("  PC%d: %.2f%%\n", k, 100 * ev[k] / total_var))
}
cat(sprintf("  Cumulative (all %d): 100.00%% (of retained)\n", M1))

# Didactic check: components are ordered by decreasing explained variance.
stopifnot(all(diff(ev) <= 1e-12))
cat("\nCheck passed: explained variance is non-increasing across PCs.\n")

# Restore X1 (undo the in-place centering/normalization).
pca1$restore()

# ==============================================================================
# Part 2: Restore Round-Trip Accuracy: Data -> PCA -> Restore -> Data
# ==============================================================================
# fit() standardizes the mapped matrix in place; restore() reverses that
# preprocessing exactly (it undoes the scaling, not the projection, so it is
# exact for any M — here M = 5 << p, as in the C++ demo).
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 2: Restore Round-Trip Accuracy\n")
cat(strrep("-", 70), "\n\n")

set.seed(123)
n2 <- 300L
p2 <- 100L
M2 <- 5L

cat(sprintf("Matrix size  : %d x %d\n", n2, p2))
cat(sprintf("Components M : %d\n\n", M2))

X2 <- matrix(rnorm(n2 * p2), nrow = n2, ncol = p2)

# Forced deep copy (`+ 0`): a plain `X2_orig <- X2` would alias the same
# buffer and be mutated together with X2 by the in-place fit().
X2_orig <- X2 + 0

pca2 <- PCA$new(X2, M = M2, center = TRUE, normalize = TRUE)
pca2$fit()

# At this point X2 holds the standardized data (in-place preprocessing).
cat("Max |X2 - X2_orig| after fit()     :",
    format(max(abs(X2 - X2_orig)), scientific = TRUE),
    " (X2 standardized in place)\n")

# restore() undoes centering and normalization on the mapped matrix.
pca2$restore()

max_error  <- max(abs(X2 - X2_orig))
mean_error <- mean(abs(X2 - X2_orig))

cat("\nReconstruction error after restore():\n")
cat("  Max Error  :", format(max_error, scientific = TRUE), "\n")
cat("  Mean Error :", format(mean_error, scientific = TRUE), "\n")

stopifnot(max_error < 1e-10)
cat("\nCheck passed: restore round-trip is exact.\n")

# ==============================================================================
# Part 3: Transform New Data: Data -> PCA -> Transform -> Scores
# ==============================================================================
# transform(X_new) applies the *training* preprocessing (means, norms) and
# projects onto the fitted loadings — new data never influences the basis.
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part 3: Transform New Data\n")
cat(strrep("-", 70), "\n\n")

set.seed(7)
n_train <- 400L
n_test  <- 50L
p3      <- 80L
M3      <- 8L

cat(sprintf("Train size   : %d x %d\n", n_train, p3))
cat(sprintf("Test size    : %d x %d\n", n_test, p3))
cat(sprintf("Components M : %d\n\n", M3))

X_train <- matrix(rnorm(n_train * p3), nrow = n_train, ncol = p3)
# Keep a pristine copy before PCA preprocesses X_train in place.
X_train_orig <- X_train + 0

X_test <- matrix(rnorm(n_test * p3), nrow = n_test, ncol = p3)

pca3 <- PCA$new(X_train, M = M3, center = TRUE, normalize = TRUE)
pca3$fit()

Z_train <- pca3$get_scores()
cat(sprintf("Train scores Z dimensions: %d x %d\n", nrow(Z_train), ncol(Z_train)))

Z_test <- pca3$transform(X_test)
cat(sprintf("Test  scores Z dimensions: %d x %d\n\n", nrow(Z_test), ncol(Z_test)))

stopifnot(nrow(Z_test) == n_test, ncol(Z_test) == M3)
cat(sprintf("Check passed: transform output shape correct (%d x %d).\n",
            n_test, M3))

# transform() on the original (unpreprocessed) training data must reproduce
# the fit scores.
Z_train_check <- pca3$transform(X_train_orig)
max_diff <- max(abs(Z_train_check - Z_train))
cat("Max diff (fit scores vs transform on same data):",
    format(max_diff, scientific = TRUE), "\n")

stopifnot(max_diff < 1e-10)
cat("Check passed: transform() reproduces the fit scores.\n")

# Restore X_train.
pca3$restore()

cat("\n", strrep("=", 70), "\n")
cat("demo_pca_01_in_memory.R complete.\n")
