# ==============================================================================
# demo_pca_01_in_memory.py
# ==============================================================================
#
# Demo illustrating PCA usage on in-memory data.
# Mirrors cpp/ml_methods/pca/demo_mlm_pca_01/demo_mlm_pca_01.cpp.
#
# Demo content:
#   - Part A: Basic fit (n = 500, p = 200, M = 10) — score matrix Z, loading
#             matrix V, and explained variance per retained component
#   - Part B: Restore round-trip (n = 300, p = 100, M = 5) — fit() preprocesses
#             X in place; restore() must return it to its original state
#   - Part C: Transform new data (n_train = 400, n_test = 50, p = 80, M = 8) —
#             project unseen observations onto the learned components and
#             verify transform() reproduces the fit scores on the training data
#
# Notes:
#   - PCA.fit(X, M) preprocesses X in place (centering / normalization), so X
#     must be a writeable, Fortran-ordered float64 array; pass a copy if the
#     original values are needed afterwards. PCA.restore() undoes the in-place
#     preprocessing of the fitted matrix.
#   - The C++ demo seeds std::mt19937; here we seed numpy's default_rng with
#     the same values. The RNG streams differ, so the printed numbers match
#     the C++ output qualitatively, not bitwise.
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods import PCA

print("\n" + "=" * 70)
print("PCA Demos (in-memory)")
print("=" * 70 + "\n")


# ==============================================================================
# Part A: Basic Fit
# ==============================================================================

print("-" * 70)
print("Part A: PCA Basic Fit")
print("-" * 70 + "\n")

n, p, M = 500, 200, 10
print(f"Matrix size  : {n} x {p}")
print(f"Components M : {M}\n")

rng = np.random.default_rng(42)
# PCA operates on the matrix in place through Eigen and therefore requires a
# writeable, Fortran-ordered (column-major) float64 array.
X = np.asfortranarray(rng.standard_normal((n, p)))

pca = PCA(center=True, normalize=True)
res = pca.fit(X, M)

# Z holds the principal-component scores: each row is an observation expressed
# in the M-dimensional component space. V holds the loadings: each column is a
# unit-length direction in the original p-dimensional variable space.
print(f"Z dimensions : {res.Z.shape[0]} x {res.Z.shape[1]}")
print(f"V dimensions : {res.V.shape[0]} x {res.V.shape[1]}")
print(f"explained_variance size: {res.explained_variance.size}\n")

# explained_variance[k] is the fraction of the total variance captured by
# component k, so the entries are decreasing and sum to at most 1 (equality
# only when M equals the full rank).
ev = res.explained_variance
assert np.all(np.diff(ev) <= 1e-12), "explained variance must be decreasing"
assert ev.sum() <= 1.0 + 1e-12, "explained variance fractions must sum to <= 1"
print(f"Total variance fraction retained by {M} components: {ev.sum():.4f}")

# Like the C++ demo, report each component as a percentage of the retained sum.
total_var = ev.sum()
print("\nExplained variance per component (% of retained):")
for k in range(M):
    print(f"  PC{k + 1}: {100.0 * ev[k] / total_var:.2f}%")
print(f"  Cumulative (all {M}): 100.00% (of retained)")

# Restore X to its original (pre-preprocessing) state.
pca.restore()


# ==============================================================================
# Part B: Restore Round-Trip Accuracy: Data -> PCA -> Restore -> Data
# ==============================================================================

print("\n" + "-" * 70)
print("Part B: PCA Restore Round-Trip Accuracy")
print("-" * 70 + "\n")

n, p, M = 300, 100, 5
print(f"Matrix size  : {n} x {p}")
print(f"Components M : {M}\n")

rng = np.random.default_rng(123)
X_orig = np.asfortranarray(rng.standard_normal((n, p)))

# Work on a copy so we can compare against the original after restore().
X = X_orig.copy(order="F")

pca = PCA()
pca.fit(X, M)

# fit() has centered/normalized X in place; restore() undoes that. The point
# is interface correctness (a round-trip on the preprocessed training matrix),
# not low-rank approximation error.
pca.restore()

max_error = np.max(np.abs(X - X_orig))
mean_error = np.mean(np.abs(X - X_orig))

print("Reconstruction error after restore():")
print(f"  Max Error  : {max_error:.6e}")
print(f"  Mean Error : {mean_error:.6e}\n")

assert max_error < 1e-10, "restore() round-trip failed"
print("Restore round-trip passed (max error < 1e-10)")


# ==============================================================================
# Part C: Transform New Data: Data -> PCA -> Transform -> Scores
# ==============================================================================

print("\n" + "-" * 70)
print("Part C: PCA Transform New Data")
print("-" * 70 + "\n")

n_train, n_test, p, M = 400, 50, 80, 8
print(f"Train size   : {n_train} x {p}")
print(f"Test size    : {n_test} x {p}")
print(f"Components M : {M}\n")

rng = np.random.default_rng(7)
X_train = np.asfortranarray(rng.standard_normal((n_train, p)))
# Keep a pristine copy before PCA preprocesses X_train in place.
X_train_orig = X_train.copy(order="F")
X_test = np.asfortranarray(rng.standard_normal((n_test, p)))

pca = PCA()
res = pca.fit(X_train, M)
print(f"Train scores Z dimensions: {res.Z.shape[0]} x {res.Z.shape[1]}")

# transform() applies the training-set centering/normalization to new data and
# projects it onto the learned loadings — no refitting takes place.
Z_test = pca.transform(X_test)
print(f"Test  scores Z dimensions: {Z_test.shape[0]} x {Z_test.shape[1]}\n")

assert Z_test.shape == (n_test, M), "transform output shape incorrect"
print(f"Transform output shape correct ({n_test} x {M})")

# Transforming the untouched training data must reproduce the fit scores.
Z_train_check = pca.transform(X_train_orig)
max_diff = np.max(np.abs(Z_train_check - res.Z))
print(f"Max diff (fit scores vs transform on same data): {max_diff:.6e}")

assert max_diff < 1e-10, "transform() does not reproduce fit scores"
print("transform() reproduces fit scores")

# Restore X_train to its original state.
pca.restore()

print("\n" + "=" * 70)
print("demo_pca_01_in_memory.py complete.")
