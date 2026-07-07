# ==============================================================================
# demo_svd_01_in_memory.py
# ==============================================================================
#
# Demo illustrating SVDSolver usage on in-memory data.
# Mirrors cpp/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01.cpp.
#
# Demo content:
#   - Part A: Direct / thin-Jacobi path (n = 300, p = 80, M = 10) — truncated
#             SVD, relative reconstruction error, top singular values, plus a
#             full-rank (M = p) sanity check
#   - Part B: Gram / kernel path (n = 100, p = 500, M = 8) — wide matrices
#             with p > 2*n route through the n x n Gram matrix
#   - Part C: Orthogonality check (n = 200, p = 150, M = 20) — U and V must
#             have orthonormal columns
#
# Notes:
#   - SVDSolver.compute(X, M) returns the rank-M truncated SVD: it keeps only
#     the M largest singular values (S), the corresponding left singular
#     vectors (U, n x M) and right singular vectors (V, p x M), so
#     U @ diag(S) @ V.T is the best rank-M approximation of X in the
#     Frobenius norm (Eckart-Young).
#   - The C++ demo seeds std::mt19937; here we seed numpy's default_rng with
#     the same values. The RNG streams differ, so the printed numbers match
#     the C++ output qualitatively, not bitwise.
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods import SVDSolver

print("\n" + "=" * 70)
print("SVDSolver Demos (in-memory)")
print("=" * 70 + "\n")


# ==============================================================================
# Part A: Direct / Thin-Jacobi Path (n >= p)
# ==============================================================================

print("-" * 70)
print("Part A: SVDSolver - Direct Path (n >= p)")
print("-" * 70 + "\n")

n, p, M = 300, 80, 10
print(f"Matrix size : {n} x {p}")
print(f"M           : {M}\n")

rng = np.random.default_rng(42)
# The solver accesses the matrix through Eigen and therefore requires a
# Fortran-ordered (column-major) float64 array.
X = np.asfortranarray(rng.standard_normal((n, p)))

solver = SVDSolver()
res = solver.compute(X, M)

print(f"U dimensions : {res.U.shape[0]} x {res.U.shape[1]}")
print(f"S size       : {res.S.size}")
print(f"V dimensions : {res.V.shape[0]} x {res.V.shape[1]}\n")

# Relative reconstruction error ||X_M - X||_F / ||X||_F. With only M = 10 of
# 80 possible components of a random Gaussian matrix, most of the energy is
# discarded, so this error is large by design — truncation keeps only the
# dominant directions.
X_approx = res.U @ np.diag(res.S) @ res.V.T
rel_error = np.linalg.norm(X_approx - X) / np.linalg.norm(X)
print(f"Relative reconstruction error ||X_M - X||_F / ||X||_F : {rel_error:.6e}")

# Singular values come back sorted in decreasing order.
print(f"\nTop {M} singular values:")
for k in range(M):
    print(f"  s{k + 1} = {res.S[k]:.4f}")
assert np.all(np.diff(res.S) <= 1e-12), "singular values must be decreasing"

# Didactic check: at full rank (M = p) the truncated SVD reproduces X exactly
# up to floating-point error.
res_full = solver.compute(X, p)
X_full = res_full.U @ np.diag(res_full.S) @ res_full.V.T
rel_error_full = np.linalg.norm(X_full - X) / np.linalg.norm(X)
print(f"\nFull-rank (M = {p}) reconstruction error: {rel_error_full:.6e}")
assert rel_error_full < 1e-10, "full-rank SVD must reconstruct X exactly"
print("Full-rank reconstruction check passed")


# ==============================================================================
# Part B: Gram / Kernel Path (p > 2*n)
# ==============================================================================

print("\n" + "-" * 70)
print("Part B: SVDSolver - Gram Path (p > 2*n)")
print("-" * 70 + "\n")

n, p, M = 100, 500, 8
print(f"Matrix size : {n} x {p}  (p > 2*n -> Gram path)")
print(f"M           : {M}\n")

rng = np.random.default_rng(7)
X = np.asfortranarray(rng.standard_normal((n, p)))

# For very wide matrices the solver works with the n x n Gram matrix X X^T
# instead of X itself — the same truncated factorization at far lower cost.
solver = SVDSolver()
res = solver.compute(X, M)

print(f"U dimensions : {res.U.shape[0]} x {res.U.shape[1]}")
print(f"S size       : {res.S.size}")
print(f"V dimensions : {res.V.shape[0]} x {res.V.shape[1]}\n")

X_approx = res.U @ np.diag(res.S) @ res.V.T
rel_error = np.linalg.norm(X_approx - X) / np.linalg.norm(X)
print(f"Relative reconstruction error ||X_M - X||_F / ||X||_F : {rel_error:.6e}")


# ==============================================================================
# Part C: Orthogonality Check
# ==============================================================================

print("\n" + "-" * 70)
print("Part C: SVDSolver - Orthogonality Check")
print("-" * 70 + "\n")

n, p, M = 200, 150, 20
print(f"Matrix size : {n} x {p}")
print(f"M           : {M}\n")

rng = np.random.default_rng(99)
X = np.asfortranarray(rng.standard_normal((n, p)))

solver = SVDSolver()
res = solver.compute(X, M)

# Singular vectors are orthonormal: U^T U = I and V^T V = I (M x M identity).
UtU_error = np.linalg.norm(res.U.T @ res.U - np.eye(M))
VtV_error = np.linalg.norm(res.V.T @ res.V - np.eye(M))

print(f"||U^T U - I||_F : {UtU_error:.6e}")
print(f"||V^T V - I||_F : {VtV_error:.6e}\n")

assert UtU_error < 1e-10 and VtV_error < 1e-10, "orthonormality check failed"
print("U and V are orthonormal")

print("\n" + "=" * 70)
print("demo_svd_01_in_memory.py complete.")
