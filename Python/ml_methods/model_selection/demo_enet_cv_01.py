# ==============================================================================
# demo_enet_cv_01.py
# ==============================================================================
#
# Demo illustrating coordinate-descent elastic-net path fitting (ElasticNet)
# and K-fold cross-validation (ElasticNetCV).
#
# Demo content:
#   Part A — ElasticNet (path fit):
#     A1. Auto-grid path for alpha = {0.0 (ridge), 0.5 (EN), 1.0 (lasso)},
#         low-dimensional (n=300, p=100, n > p); each auto path is also
#         reproduced via fit_grid() at the same explicit lambda grid.
#     A2. High-dimensional path (n=200, p=500, p > n) for alpha = 0.0 and 1.0.
#   Part B — ElasticNetCV (K-fold CV):
#     B1. CV lambda selection (lambda.min / lambda.1se) for each alpha,
#         low-dimensional (n=300, p=100, 10-fold).
#     B2. CV for pure ridge (alpha=0) as used by TRexGVS; shows the
#         lambda_2_lars = lambda.1se * p / 2 conversion.
#
# Mirrors: cpp/ml_methods/model_selection/demo_mlm_ms_02_enet_cv_ccd/
#          demo_mlm_ms_02_enet_cv_ccd.cpp
#
# Notes on the port:
#   - Same n / p / alpha / fold values and data seeds as the C++ demo;
#     NumPy's PCG64 stream differs from std::mt19937, so results match
#     qualitatively, not bitwise.
#   - ElasticNet replicates glmnet's pathwise CCD with glmnet's fdev/devmax
#     early termination, so the auto grid may stop before n_lambda points
#     for alpha > 0. ElasticNetCV builds the full-data lambda grid once and
#     reuses it across folds (glmnet semantics); its grid is stored in
#     glmnet's descending order (largest penalty first), unlike RidgeCV's
#     ascending grid.
#   - Arrays must be float64 and Fortran-ordered (column-major);
#     np.asfortranarray() provides the required layout.
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods import ElasticNet, ElasticNetCV

# ==============================================================================

print("\n" + "=" * 70)
print("ElasticNet / ElasticNetCV Demo Suite (CCD path fit + K-fold CV)")
print("=" * 70 + "\n")

# ==============================================================================
# Shared helpers
# ==============================================================================


def generate_data(n, p, sigma, seed, n_active=10):
    """Synthetic regression data with known sparse signal: n_active coefficients
    equal to 1.0 at evenly spaced column indices, y = X beta_true + N(0, sigma^2).
    Mirrors the C++ generate_data() helper."""
    rng = np.random.default_rng(seed)
    X = np.asfortranarray(rng.standard_normal((n, p)))
    beta_true = np.zeros(p)
    stride = max(1, p // n_active)
    beta_true[np.arange(n_active) * stride] = 1.0
    y = X @ beta_true + sigma * rng.standard_normal(n)
    return X, y, beta_true


def print_path_summary(en, step=10):
    """Condensed coefficient path from an ElasticNet fit.

    Shows every `step`-th lambda along with the number of non-zero
    coefficients (l0-norm), the L1-norm of the coefficient vector, and the
    deviance ratio (fraction of null deviance explained — glmnet's dev.ratio,
    which increases as lambda decreases and the fit becomes less penalized)."""
    lam = np.asarray(en.lambdas())
    B   = np.asarray(en.coef())          # p x K coefficient path
    dev = np.asarray(en.dev_ratio())
    K   = len(lam)

    tail = "" if en.converged() else f" [{en.n_nonconverged()} non-converged]"
    print(f"  Path: {K} lambdas fitted{tail}")
    print(f"  {'lambda':<12}{'nnz':<8}{'L1(beta)':<12}{'dev.ratio':<10}")

    def print_row(k):
        nnz = int(np.sum(np.abs(B[:, k]) > 1e-9))
        print(f"  {lam[k]:<12.5f}{nnz:<8d}"
              f"{np.sum(np.abs(B[:, k])):<12.5f}{dev[k]:<10.5f}")

    for k in range(0, K, step):
        print_row(k)
    # Always print the last point
    if (K - 1) % step != 0:
        print_row(K - 1)
    print()


def print_cv_summary(cv, step=10):
    """Condensed CV summary from an ElasticNetCV fit:
      - lambda.min: grid value minimizing the CV mean-squared error
      - lambda.1se: largest lambda whose CV-MSE is within one standard error
        of the minimum (more parsimonious / more heavily regularized choice)"""
    lam  = np.asarray(cv.lambdas())
    mse  = np.asarray(cv.cv_mse())
    sem  = np.asarray(cv.cv_sem())
    K    = len(lam)
    imin = cv.index_min()
    i1se = cv.index_1se()

    print(f"  lambda.min : {cv.cv_min():.5f}  "
          f"(CV-MSE={mse[imin]:.5f} +/- {sem[imin]:.5f})")
    print(f"  lambda.1se : {cv.cv_1se():.5f}  "
          f"(CV-MSE={mse[i1se]:.5f} +/- {sem[i1se]:.5f})\n")

    print(f"  CV-MSE curve (every {step} steps, lambda descending):")
    print(f"    {'lambda':<14}{'CV-MSE':<12}{'SEM':<12}marker")

    def print_row(k, marker):
        print(f"    {lam[k]:<14.5f}{mse[k]:<12.5f}{sem[k]:<12.5f}{marker}")

    for k in range(0, K, step):
        marker = "  <-- min" if k == imin else "  <-- 1se" if k == i1se else ""
        print_row(k, marker)
    if imin % step != 0:
        print_row(imin, "  <-- min")
    if i1se % step != 0 and i1se != imin:
        print_row(i1se, "  <-- 1se")
    print()


def fit_path_suite(title, n, p, seed, alphas, sigma=1.0):
    """Fits an auto-grid path per alpha, prints the condensed summary, and
    applies two didactic checks:
      - the deviance ratio is non-decreasing along the (descending-lambda)
        path — less penalty can only improve the training fit;
      - fit_grid() at the auto-generated grid reproduces the auto path
        exactly (fit() == fit_grid() when handed the same lambdas)."""
    print("-" * 70)
    print(title)
    print("-" * 70 + "\n")

    X, y, beta_true = generate_data(n, p, sigma, seed)
    print(f"Data: n={n}, p={p}, active coefs={int(np.sum(beta_true != 0))}\n")

    for alpha in alphas:
        # alpha mixes the penalties: alpha=1 is the pure L1 lasso (sparse
        # paths), alpha=0 pure L2 ridge (dense shrinkage), 0 < alpha < 1
        # the elastic net in between.
        print(f"--- alpha={alpha:g} ---")
        en = ElasticNet().fit(X, y, alpha)
        print_path_summary(en, step=10)

        dev = np.asarray(en.dev_ratio())
        assert np.all(np.diff(dev) >= -1e-12), \
            "dev.ratio must be non-decreasing along the path"

        # Explicit-grid solver check: refitting at the auto grid via
        # fit_grid() must reproduce the same coefficient path.
        lam_grid = np.asarray(en.lambdas(), dtype=np.float64)
        en_grid = ElasticNet().fit_grid(X, y, lam_grid, alpha)
        max_diff = np.max(np.abs(np.asarray(en.coef())
                                 - np.asarray(en_grid.coef())))
        assert max_diff < 1e-12
        print(f"  fit_grid() reproduction check "
              f"(max |coef diff| = {max_diff:.2e}): PASS\n")


# ==============================================================================
# Part A1: ElasticNet path — low-dimensional (n=300, p=100, n > p)
# ==============================================================================

fit_path_suite("Part A1: ElasticNet path | Low-Dimensional (n=300, p=100)",
               n=300, p=100, seed=101, alphas=(0.0, 0.5, 1.0))

# ==============================================================================
# Part A2: ElasticNet path — high-dimensional (n=200, p=500, p > n)
# ==============================================================================

fit_path_suite("Part A2: ElasticNet path | High-Dimensional (n=200, p=500)",
               n=200, p=500, seed=202, alphas=(0.0, 1.0))

# ==============================================================================
# Part B1: ElasticNetCV — low-dimensional CV (n=300, p=100, 10-fold)
# ==============================================================================

print("-" * 70)
print("Part B1: ElasticNetCV | Low-Dimensional CV (n=300, p=100, 10-fold)")
print("-" * 70 + "\n")

n, p, sigma, n_folds = 300, 100, 1.0, 10
X, y, beta_true = generate_data(n, p, sigma, seed=303)
print(f"Data: n={n}, p={p}, folds={n_folds}, "
      f"active coefs={int(np.sum(beta_true != 0))}\n")

for alpha in (0.0, 0.5, 1.0):
    print(f"--- alpha={alpha:g} ---")
    cv = ElasticNetCV().fit(X, y, alpha=alpha, n_folds=n_folds, seed=0)
    print_cv_summary(cv, step=10)

    # The CV curve is evaluated at every grid point (the grid may hold fewer
    # than n_lambda points when glmnet's early termination kicks in).
    assert len(cv.cv_mse()) == len(cv.lambdas())
    assert len(cv.cv_sem()) == len(cv.lambdas())
    assert cv.cv_1se() >= cv.cv_min()
    print("  Ordering check (lambda.min <= lambda.1se): PASS\n")

# ==============================================================================
# Part B2: ElasticNetCV — ridge CV + lambda_2_lars conversion (n=300, p=200)
# ==============================================================================
# TRexGVS converts the glmnet-scale ridge lambda via:
#   lambda_2_lars = lambda.1se * p / 2
# mirroring R's `lm_dummy` and `cv.glmnet(alpha=0)`.
# ==============================================================================

print("-" * 70)
print("Part B2: ElasticNetCV | Ridge CV + lambda_2_lars Conversion")
print("-" * 70 + "\n")

n, p, sigma, n_folds = 300, 200, 1.0, 10
X, y, beta_true = generate_data(n, p, sigma, seed=404)
print(f"Data: n={n}, p={p}, folds={n_folds}, "
      f"active coefs={int(np.sum(beta_true != 0))}\n")

cv = ElasticNetCV().fit(X, y, alpha=0.0, n_folds=n_folds, seed=0)
print_cv_summary(cv, step=10)

assert cv.cv_1se() >= cv.cv_min()
lambda2_min = cv.cv_min() * p / 2.0
lambda2_1se = cv.cv_1se() * p / 2.0

print(f"  lambda_2_lars (min) = cv_min  * p/2 = "
      f"{cv.cv_min():.6f} * {p}/2 = {lambda2_min:.6f}")
print(f"  lambda_2_lars (1se) = cv_1se  * p/2 = "
      f"{cv.cv_1se():.6f} * {p}/2 = {lambda2_1se:.6f}\n")
print("  (This is the value TRexGVS uses internally via CV_1SE_CCD.)\n")

print("=" * 70)
print("demo_enet_cv_01.py complete.")
