# ==============================================================================
# demo_standardization_01_in_memory.R
# ==============================================================================
#
# Demo illustrating standardization methods on in-memory data.
#
# Demo content:
#   - Z-scaling via ZScoreScaler (center + scale to unit variance)
#   - L2-norm scaling via LpNormScaler (center + scale to unit L2 norm)
#   - Forward and inverse transformations
#
# ==============================================================================

library(TRexSelector)

# ==============================================================================

set.seed(42)
n <- 8
p <- 4

X <- matrix(rnorm(n * p, mean = 3, sd = 2), nrow = n, ncol = p)

summarize_matrix <- function(mat) {
  data.frame(
    mean = round(colMeans(mat), 6),
    sd = round(apply(mat, 2L, sd), 6),
    l2_norm = round(apply(mat, 2L, function(x) sqrt(sum(x^2))), 6)
  )
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Standardization Demos (in-memory)\n")
cat(strrep("=", 70), "\n\n")

cat("Original in-memory matrix X:\n")
print(round(X, 4))
cat("\nColumn summary of original X:\n")
print(summarize_matrix(X))

# ==============================================================================
# Part A: Z-scaling
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part A: Z-Scaling (ZScoreScaler)\n")
cat(strrep("-", 70), "\n\n")

z_scaler <- ZScoreScaler$new(with_mean = TRUE, with_std = TRUE)
z_scaler$fit(X)
X_z <- z_scaler$transform_inplace(X + 0)

cat("Z-scaled matrix:\n")
print(round(X_z, 4))
cat("\nColumn summary after z-scaling:\n")
print(summarize_matrix(X_z))

X_z_back <- z_scaler$inverse_transform_inplace(X_z + 0)
cat("\nMax reconstruction error after inverse z-scaling:",
    format(max(abs(X - X_z_back)), scientific = TRUE), "\n")

# ==============================================================================
# Part B: L2-norm scaling
# ==============================================================================

cat("\n", strrep("-", 70), "\n")
cat("Part B: L2-Norm Scaling (LpNormScaler, norm_type = 2)\n")
cat(strrep("-", 70), "\n\n")

l2_scaler <- LpNormScaler$new(norm_type = 2L, with_mean = TRUE)
l2_scaler$fit(X)
X_l2 <- l2_scaler$transform_inplace(X + 0)

cat("L2-scaled matrix:\n")
print(round(X_l2, 4))
cat("\nColumn summary after L2 scaling:\n")
print(summarize_matrix(X_l2))

X_l2_back <- l2_scaler$inverse_transform_inplace(X_l2 + 0)
cat("\nMax reconstruction error after inverse L2 scaling:",
    format(max(abs(X - X_l2_back)), scientific = TRUE), "\n")

cat("\n", strrep("=", 70), "\n")
cat("demo_standardization_01_in_memory.R complete.\n")
