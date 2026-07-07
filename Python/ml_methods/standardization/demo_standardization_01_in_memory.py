# ==============================================================================
# demo_standardization_01_in_memory.py
# ==============================================================================
#
# Demo illustrating standardization methods on in-memory data.
# Mirrors R/ml_methods/standardization/demo_standardization_01_in_memory.R.
#
# Demo content:
#   - Z-scaling via ZScoreScaler (center + scale to unit variance)
#   - L2-norm scaling via LpNormScaler (center + scale to unit L2 norm)
#   - Forward and inverse transformations
#
# Notes:
#   transform_inplace() mutates the supplied array in place and returns None.
#   We therefore copy first (order="F" keeps the required column-major layout,
#   mirroring R's X + 0) and pass the copy, preserving the original X for the
#   reconstruction check.
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods import LpNormScaler, NormType, ZScoreScaler

# ==============================================================================
# Global Parameters and Data Generation
# ==============================================================================

rng = np.random.default_rng(42)

n = 8
p = 4
# The scalers operate on the matrix in-place through Eigen and therefore
# require a writeable, Fortran-ordered (column-major) float64 array.
X = np.asfortranarray(rng.normal(3.0, 2.0, (n, p)))


def summarize_matrix(mat):
    """Print column-wise mean, std (ddof=1), and L2-norm — mirrors R summarize_matrix()."""
    means   = np.mean(mat, axis=0)
    stds    = np.std(mat, axis=0, ddof=1)
    l2norms = np.sqrt(np.sum(mat ** 2, axis=0))
    header  = f"  {'col':>4}  {'mean':>10}  {'std':>10}  {'l2_norm':>10}"
    print(header)
    for j in range(mat.shape[1]):
        print(f"  {j:>4}  {means[j]:>10.6f}  {stds[j]:>10.6f}  {l2norms[j]:>10.6f}")


print("\n" + "=" * 70)
print("Standardization Demos (in-memory)")
print("=" * 70 + "\n")

print("Original in-memory matrix X:")
print(np.round(X, 4))
print("\nColumn summary of original X:")
summarize_matrix(X)


# ==============================================================================
# Part A: Z-scaling (ZScoreScaler)
# ==============================================================================

print("\n" + "-" * 70)
print("Part A: Z-Scaling (ZScoreScaler)")
print("-" * 70 + "\n")

z_scaler = ZScoreScaler(center=True, scale=True)
z_scaler.fit(X)

# Copy first (mirrors R's X + 0), then transform the copy in place.
X_z = X.copy(order="F")
z_scaler.transform_inplace(X_z)

print("Z-scaled matrix:")
print(np.round(X_z, 4))
print("\nColumn summary after z-scaling:")
summarize_matrix(X_z)

X_z_back = X_z.copy(order="F")
z_scaler.inverse_transform_inplace(X_z_back)
recon_err = np.max(np.abs(X - X_z_back))
print(f"\nMax reconstruction error after inverse z-scaling: {recon_err:.2e}")


# ==============================================================================
# Part B: L2-norm scaling (LpNormScaler)
# ==============================================================================

print("\n" + "-" * 70)
print("Part B: L2-Norm Scaling (LpNormScaler, norm_type = L2)")
print("-" * 70 + "\n")

l2_scaler = LpNormScaler(NormType.L2, True)
l2_scaler.fit(X)

X_l2 = X.copy(order="F")
l2_scaler.transform_inplace(X_l2)

print("L2-scaled matrix:")
print(np.round(X_l2, 4))
print("\nColumn summary after L2 scaling:")
summarize_matrix(X_l2)

X_l2_back = X_l2.copy(order="F")
l2_scaler.inverse_transform_inplace(X_l2_back)
recon_err_l2 = np.max(np.abs(X - X_l2_back))
print(f"\nMax reconstruction error after inverse L2 scaling: {recon_err_l2:.2e}")

print("\n" + "=" * 70)
print("demo_standardization_01_in_memory.py complete.")
