# ==============================================================================
# demo_mlm_svd_01.py
# ==============================================================================
#
# Demonstration of SVDSolver usage.
#
# Mirrors cpp/ml_methods/svd/demo_mlm_svd_01/demo_mlm_svd_01.cpp.
#
# Three demos which illustrate the truncated SVD solver:
#
#   Demo 1 | Direct path   : n = 300, p = 80, M = 10 — thin-Jacobi SVD for
#                            n >= p, relative reconstruction error, top
#                            singular values (plus a full-rank M = p sanity
#                            check, a Python-side addition).
#
#   Demo 2 | Gram path     : n = 100, p = 500, M = 8 — wide matrices with
#                            p > 2*n route through the n x n Gram matrix.
#
#   Demo 3 | Orthogonality : n = 200, p = 150, M = 20 — U and V must have
#                            orthonormal columns.
#
# Notes:
#   - SVDSolver.compute(X, M) returns the rank-M truncated SVD: it keeps only
#     the M largest singular values (S), the corresponding left singular
#     vectors (U, n x M) and right singular vectors (V, p x M), so
#     U @ diag(S) @ V.T is the best rank-M approximation of X in the
#     Frobenius norm (Eckart-Young).
#   - The C++ demo seeds std::mt19937; here numpy's default_rng is seeded with
#     the same values. The RNG streams differ, so the printed numbers match
#     the C++ output qualitatively, not bitwise.
#
# ==============================================================================

import numpy as np

from trex_selector_neo.ml_methods import SVDSolver


def print_section_header(title):
    print("=" * 78)
    print(f" {title}")
    print("=" * 78 + "\n")


# ==============================================================================
# Demo 1: Direct / Thin-Jacobi Path (n >= p)
# ==============================================================================

def demo_svd_direct():
    print_section_header("SVDSolver - Direct Path (n >= p)")

    # Step 1. Generate random data
    n, p, M = 300, 80, 10
    print(f"Matrix size : {n} x {p}")
    print(f"M           : {M}\n")

    rng = np.random.default_rng(42)
    # The solver accesses the matrix through Eigen and therefore requires a
    # Fortran-ordered (column-major) float64 array.
    X = np.asfortranarray(rng.standard_normal((n, p)))

    # Step 2. Compute SVD
    solver = SVDSolver()
    res = solver.compute(X, M)

    # Step 3. Print dimensions
    print(f"U dimensions : {res.U.shape[0]} x {res.U.shape[1]}")
    print(f"S size       : {res.S.size}")
    print(f"V dimensions : {res.V.shape[0]} x {res.V.shape[1]}\n")

    # Step 4. Reconstruction error: ||X_M - X||_F / ||X||_F. With only M = 10
    #         of 80 possible components of a random Gaussian matrix, most of
    #         the energy is discarded, so this error is large by design —
    #         truncation keeps only the dominant directions.
    X_approx = res.U @ np.diag(res.S) @ res.V.T
    rel_error = np.linalg.norm(X_approx - X) / np.linalg.norm(X)
    print(f"Relative reconstruction error ||X_M - X||_F / ||X||_F : {rel_error:.6e}")

    # Step 5. Singular values (returned sorted in decreasing order)
    print(f"\nTop {M} singular values:")
    for k in range(M):
        print(f"  s{k + 1} = {res.S[k]:.4f}")

    # Python-side sanity check (not in the cpp demo): at full rank (M = p) the
    # truncated SVD reproduces X exactly up to floating-point error.
    res_full = solver.compute(X, p)
    X_full = res_full.U @ np.diag(res_full.S) @ res_full.V.T
    rel_error_full = np.linalg.norm(X_full - X) / np.linalg.norm(X)
    print(f"\nFull-rank (M = {p}) reconstruction error: {rel_error_full:.6e}")
    assert rel_error_full < 1e-10, "full-rank SVD must reconstruct X exactly"
    print("Full-rank reconstruction check passed")
    print("\n")


# ==============================================================================
# Demo 2: Gram / Kernel Path (p > 2*n)
# ==============================================================================

def demo_svd_gram():
    print_section_header("SVDSolver - Gram Path (p > 2*n)")

    # Step 1. Generate wide matrix
    n, p, M = 100, 500, 8
    print(f"Matrix size : {n} x {p}  (p > 2*n -> Gram path)")
    print(f"M           : {M}\n")

    rng = np.random.default_rng(7)
    X = np.asfortranarray(rng.standard_normal((n, p)))

    # Step 2. Compute SVD — for very wide matrices the solver works with the
    #         n x n Gram matrix X X^T instead of X itself: the same truncated
    #         factorization at far lower cost.
    solver = SVDSolver()
    res = solver.compute(X, M)

    # Step 3. Print dimensions
    print(f"U dimensions : {res.U.shape[0]} x {res.U.shape[1]}")
    print(f"S size       : {res.S.size}")
    print(f"V dimensions : {res.V.shape[0]} x {res.V.shape[1]}\n")

    # Step 4. Reconstruction error
    X_approx = res.U @ np.diag(res.S) @ res.V.T
    rel_error = np.linalg.norm(X_approx - X) / np.linalg.norm(X)
    print(f"Relative reconstruction error ||X_M - X||_F / ||X||_F : {rel_error:.6e}")
    print("\n")


# ==============================================================================
# Demo 3: Orthogonality Check
# ==============================================================================

def demo_svd_orthogonality():
    print_section_header("SVDSolver - Orthogonality Check")

    n, p, M = 200, 150, 20
    print(f"Matrix size : {n} x {p}")
    print(f"M           : {M}\n")

    rng = np.random.default_rng(99)
    X = np.asfortranarray(rng.standard_normal((n, p)))

    solver = SVDSolver()
    res = solver.compute(X, M)

    # Singular vectors are orthonormal: U^T U = I and V^T V = I (M x M).
    UtU_error = np.linalg.norm(res.U.T @ res.U - np.eye(M))
    VtV_error = np.linalg.norm(res.V.T @ res.V - np.eye(M))

    print(f"||U^T U - I||_F : {UtU_error:.6e}")
    print(f"||V^T V - I||_F : {VtV_error:.6e}\n")

    if UtU_error < 1e-10 and VtV_error < 1e-10:
        print("✓ U and V are orthonormal")
    else:
        print("✗ Orthonormality check FAILED")
    print("\n")


# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    demo_svd_direct()
    demo_svd_gram()
    demo_svd_orthogonality()
