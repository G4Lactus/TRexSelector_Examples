# ==============================================================================
# demo_standardization_02_mmap.py
# ==============================================================================
#
# Demo illustrating standardization methods on memory-mapped data.
# Mirrors R/ml_methods/standardization/demo_standardization_02_mmap.R.
#
# Demo content:
#   - Converting in-memory matrix to a memory-mapped file
#   - Z-scaling via ZScoreScaler on mmap-backed data
#   - L2-norm scaling via LpNormScaler on mmap-backed data
#   - Forward and inverse transformations
#
# Key note:
#   X_mmap.to_numpy() returns a zero-copy view into the mmap buffer.
#   Calling transform_inplace() on that view directly would mutate the binary
#   file on disk.  Instead we take X_loaded = X_mmap.to_numpy().copy(order="F"),
#   which mirrors R's X_loaded <- as.matrix(X_mmap): both produce an independent
#   in-memory (column-major) copy that the scaler is free to modify in place.
#   transform_inplace() mutates its argument and returns None, so we always
#   copy first and pass the copy.
#
# ==============================================================================

import os
import tempfile

import numpy as np
from trex_selector_neo.ml_methods import LpNormScaler, NormType, ZScoreScaler
from trex_selector_neo.utils import numpy_to_memmap

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


# Create temporary mmap backing file
fd, mmap_path = tempfile.mkstemp(suffix=".bin")
os.close(fd)

try:
    # Write X to the mmap file; obtain a zero-copy view; then take a copy.
    # The copy mirrors R's `as.matrix(X_mmap)` and is required so that
    # transform_inplace() does not mutate the binary file on disk.
    X_mmap   = numpy_to_memmap(mmap_path, X)
    X_loaded = X_mmap.to_numpy().copy(order="F")

    print("\n" + "=" * 70)
    print("Standardization Demos (memory-mapped)")
    print("=" * 70 + "\n")

    print("Original memory-mapped matrix X:")
    print(np.round(X_loaded, 4))
    print("\nColumn summary of original mmap-backed X:")
    summarize_matrix(X_loaded)

    # ==========================================================================
    # Part A: Z-scaling (ZScoreScaler)
    # ==========================================================================

    print("\n" + "-" * 70)
    print("Part A: Z-Scaling (ZScoreScaler) on mmap-backed data")
    print("-" * 70 + "\n")

    z_scaler = ZScoreScaler(center=True, scale=True)
    z_scaler.fit(X_loaded)

    X_z = X_loaded.copy(order="F")
    z_scaler.transform_inplace(X_z)

    print("Z-scaled matrix:")
    print(np.round(X_z, 4))
    print("\nColumn summary after z-scaling:")
    summarize_matrix(X_z)

    X_z_back = X_z.copy(order="F")
    z_scaler.inverse_transform_inplace(X_z_back)
    recon_err = np.max(np.abs(X_loaded - X_z_back))
    print(f"\nMax reconstruction error after inverse z-scaling: {recon_err:.2e}")

    # ==========================================================================
    # Part B: L2-norm scaling (LpNormScaler)
    # ==========================================================================

    print("\n" + "-" * 70)
    print("Part B: L2-Norm Scaling (LpNormScaler, norm_type = L2) on mmap-backed data")
    print("-" * 70 + "\n")

    l2_scaler = LpNormScaler(NormType.L2, True)
    l2_scaler.fit(X_loaded)

    X_l2 = X_loaded.copy(order="F")
    l2_scaler.transform_inplace(X_l2)

    print("L2-scaled matrix:")
    print(np.round(X_l2, 4))
    print("\nColumn summary after L2 scaling:")
    summarize_matrix(X_l2)

    X_l2_back = X_l2.copy(order="F")
    l2_scaler.inverse_transform_inplace(X_l2_back)
    recon_err_l2 = np.max(np.abs(X_loaded - X_l2_back))
    print(f"\nMax reconstruction error after inverse L2 scaling: {recon_err_l2:.2e}")

finally:
    try:
        os.unlink(mmap_path)
    except OSError:
        pass

print("\n" + "=" * 70)
print("demo_standardization_02_mmap.py complete.")
