# ==============================================================================
# demo_ridge_cv_01.py
# ==============================================================================
#
# Demo illustrating ridge-regression K-fold cross-validation via RidgeCV
# (SVD-based solver with per-fold centering and column L2-normalization).
# Mirrors R/ml_methods/model_selection/demo_ridge_cv_01.R.
#
# Demo content:
#   - Scenario 1: Low-dimensional  (n=1000, p=500,  n > p,  10-fold)
#   - Scenario 2: High-dimensional (n=300,  p=500,  p > n,  10-fold)
#   - Scenario 3: Very high-dim    (n=300,  p=1000, p >> n,  5-fold)
#   - Each scenario reports the lambda grid endpoints, the CV-MSE curve,
#     and the lambda.min / lambda.1se selections.
#
# Mirrors: cpp/ml_methods/model_selection/demo_mlm_ms_01_ridge_cv_svd/
#          demo_mlm_ms_01_ridge_cv_svd.cpp
#
# Notes on the port:
#   - Same n / p / fold values and data seeds as the C++ and R demos; NumPy's
#     PCG64 stream differs from std::mt19937 and from R's Mersenne-Twister,
#     so results match qualitatively, not bitwise.
#   - The C++ demo has a 4th, memory-mapped scenario (the design matrix is
#     read from disk via a zero-copy Eigen map). RidgeCV's Python binding
#     only accepts in-memory arrays, so that scenario is omitted here.
#   - RidgeCV requires float64, Fortran-ordered (column-major) arrays;
#     np.asfortranarray() provides the required layout.
#
# ==============================================================================

import numpy as np
from trex_selector_neo.ml_methods import RidgeCV

# ==============================================================================

print("\n" + "=" * 70)
print("RidgeCV Demo Suite (K-fold CV for the ridge penalty)")
print("=" * 70 + "\n")

# ==============================================================================
# Shared helpers
# ==============================================================================


def generate_regression_data(n, p, sigma, seed, n_active=10):
    """Synthetic regression data with known sparse signal: n_active coefficients
    equal to 1.0 at evenly spaced column indices, y = X beta_true + N(0, sigma^2).
    Mirrors the C++ generate_regression_data() helper."""
    rng = np.random.default_rng(seed)
    X = np.asfortranarray(rng.standard_normal((n, p)))
    beta_true = np.zeros(p)
    stride = max(1, p // n_active)
    beta_true[np.arange(n_active) * stride] = 1.0
    y = X @ beta_true + sigma * rng.standard_normal(n)
    return X, y, beta_true


def print_cv_summary(cv, step_size=10):
    """Condensed summary of a fitted RidgeCV object:
      - lambda.min: grid value minimizing the CV mean-squared error
      - lambda.1se: largest lambda whose CV-MSE is within one standard error
        of the minimum (more parsimonious / more heavily regularized choice)
    The lambda grid is stored in ascending order; the curve is printed every
    `step_size` grid points with the min/1se positions marked."""
    lambdas = np.asarray(cv.lambdas())
    mse     = np.asarray(cv.cv_mse())
    sem     = np.asarray(cv.cv_sem())
    K       = len(lambdas)
    i_min   = cv.index_min()
    i_1se   = cv.index_1se()

    print(f"  lambda.min : {cv.cv_min():.4f}  "
          f"(CV-MSE={mse[i_min]:.4f} +/- {sem[i_min]:.4f})")
    print(f"  lambda.1se : {cv.cv_1se():.4f}  "
          f"(CV-MSE={mse[i_1se]:.4f} +/- {sem[i_1se]:.4f})\n")

    print(f"  CV-MSE curve (every {step_size} steps, lambda ascending):")
    print(f"    {'lambda':<14}{'CV-MSE':<12}{'SEM':<12}marker")

    def print_row(k, marker):
        print(f"    {lambdas[k]:<14.4f}{mse[k]:<12.4f}{sem[k]:<12.4f}{marker}")

    shown = range(0, K, step_size)
    for k in shown:
        marker = "  <-- min" if k == i_min else "  <-- 1se" if k == i_1se else ""
        print_row(k, marker)
    if i_min not in shown:
        print_row(i_min, "  <-- min")
    if i_1se not in shown and i_1se != i_min:
        print_row(i_1se, "  <-- 1se")
    print()


def run_cv_scenario(title, n, p, n_folds, seed, sigma=1.0):
    """Runs one CV scenario end to end and applies the ordering check
    lambda.1se >= lambda.min (the 1se rule always regularizes at least as
    strongly as the pure minimizer)."""
    print("-" * 70)
    print(title)
    print("-" * 70 + "\n")

    X, y, beta_true = generate_regression_data(n, p, sigma, seed)
    print(f"Data: n={n}, p={p}, folds={n_folds}, "
          f"active coefs={int(np.sum(beta_true != 0))}\n")

    # Defaults mirror the C++ fit(): 100-lambda grid anchored at the largest
    # useful penalty, spanning a 1000x range; `seed` fixes fold assignment.
    cv = RidgeCV()
    cv.fit(X, y, n_folds=n_folds, n_lambda=100, lambda_ratio=1000.0, seed=seed)

    print_cv_summary(cv, step_size=10)

    assert cv.cv_1se() >= cv.cv_min()
    print("  Ordering check (lambda.min <= lambda.1se): PASS\n")
    return cv


# ==============================================================================
# Scenario 1: Low-dimensional (n=1000, p=500, n > p, 10-fold)
# ==============================================================================

run_cv_scenario("Scenario 1: Low-Dimensional (n=1000, p=500, 10-fold)",
                n=1000, p=500, n_folds=10, seed=111)

# ==============================================================================
# Scenario 2: High-dimensional (n=300, p=500, p > n, 10-fold)
# ==============================================================================

run_cv_scenario("Scenario 2: High-Dimensional (n=300, p=500, 10-fold)",
                n=300, p=500, n_folds=10, seed=222)

# ==============================================================================
# Scenario 3: Very high-dimensional (n=300, p=1000, p >> n, 5-fold)
# ==============================================================================
# 5 folds keep each fold's training set large enough relative to n.
# ==============================================================================

run_cv_scenario("Scenario 3: Very High-Dimensional (n=300, p=1000, 5-fold)",
                n=300, p=1000, n_folds=5, seed=333)

print("=" * 70)
print("demo_ridge_cv_01.py complete.")
