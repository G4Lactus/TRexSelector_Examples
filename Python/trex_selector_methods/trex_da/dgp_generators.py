# ==============================================================================
# dgp_generators.py
# ==============================================================================
#
# Data-generating processes for the DA-TRex demos, with SNR-calibrated noise
# (SNR = Var(signal) / Var(noise), sample variance with the n-1 denominator to
# match the C++ make_snr_response and the R ports).
#
# Support indices are 0-based throughout (Python convention).
#
# Mirrors cpp/trex_selector_methods/trex_da/dgp_generators.hpp and the R
# dgp_*_snr.R files.
# ==============================================================================

import numpy as np


def _snr_response(X, beta, snr, rng):
    """y = X @ beta + noise, with noise_sigma = sqrt(var(signal) / snr)."""
    signal = X @ beta
    noise_sigma = np.sqrt(np.var(signal, ddof=1) / snr)
    y = signal + rng.standard_normal(X.shape[0]) * noise_sigma
    return y


def dgp_ar1_snr(n, p, support, amplitude=3.0, rho=0.7, snr=2.0, seed=None):
    """AR(1) Toeplitz DGP: Sigma[j, k] = rho^|j-k|, SNR-calibrated noise.

    X[:, 0] = eta_0;  X[:, j] = rho * X[:, j-1] + sqrt(1 - rho^2) * eta_j.
    """
    assert abs(rho) < 1 and snr > 0
    rng = np.random.default_rng(seed)
    scale = np.sqrt(1.0 - rho * rho)

    X = np.empty((n, p))
    X[:, 0] = rng.standard_normal(n)
    for j in range(1, p):
        X[:, j] = rho * X[:, j - 1] + scale * rng.standard_normal(n)

    beta = np.zeros(p)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), snr=snr, rho=rho)


def dgp_bt_snr(n, p, support, amplitude=3.0, n_blocks=10, rho_within=0.8,
               rho_between=0.1, snr=2.0, seed=None):
    """Two-level factor-model block covariance: Sigma = rho_within within a
    block, rho_between across blocks. With rho_within == rho_between this is
    pure equicorrelation (compound symmetry). SNR-calibrated noise.

    X = sqrt(rho_between)*f0 + sqrt(rho_within-rho_between)*f_block + sqrt(1-rho_within)*eta
    """
    assert p % n_blocks == 0
    assert 0.0 <= rho_between <= rho_within < 1.0 and snr > 0.0

    rng = np.random.default_rng(seed)
    block_size = p // n_blocks
    labels = np.repeat(np.arange(n_blocks), block_size)   # length p

    a = np.sqrt(rho_between)
    b = np.sqrt(rho_within - rho_between)
    c = np.sqrt(1.0 - rho_within)

    f0  = rng.standard_normal(n)
    fB  = rng.standard_normal((n, n_blocks))
    eta = rng.standard_normal((n, p))
    X = a * f0[:, None] + b * fB[:, labels] + c * eta

    beta = np.zeros(p)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), n_blocks=n_blocks,
                rho_within=rho_within, rho_between=rho_between, snr=snr)


def dgp_nn_snr(n, p, support, amplitude=3.0, kappa=3, rho=0.7, snr=2.0, seed=None):
    """Nearest-Neighbours / MA(kappa) banded DGP: each column is a causal MA(kappa)
    convolution of a shared innovation sheet, giving Sigma[j,k] != 0 only for
    |j-k| <= kappa. theta[l] = rho^l, normalised to unit variance.
    """
    assert kappa >= 1 and 0 <= rho < 1 and snr > 0
    rng = np.random.default_rng(seed)

    theta = rho ** np.arange(kappa + 1)
    theta = theta / np.sqrt(np.sum(theta ** 2))
    theta_rev = theta[::-1]

    eta = rng.standard_normal((n, p + kappa))       # padded innovation sheet
    X = np.empty((n, p))
    for j in range(p):
        X[:, j] = eta[:, j:j + kappa + 1] @ theta_rev

    beta = np.zeros(p)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), kappa=kappa, rho=rho, snr=snr)


def dgp_ar1_block_snr(n, p, support, amplitude=3.0, n_blocks=5, rho=0.7, snr=2.0, seed=None):
    """Block-diagonal AR(1) DGP: p predictors split into n_blocks blocks; within
    each block the columns follow an AR(1) recursion (Toeplitz rho^|j-k|), zero
    correlation between blocks. SNR-calibrated noise.
    """
    assert p % n_blocks == 0 and abs(rho) < 1 and snr > 0
    rng = np.random.default_rng(seed)
    block_size = p // n_blocks
    innov_sd = np.sqrt(1.0 - rho * rho)

    X = np.zeros((n, p))
    for m in range(n_blocks):
        cs = m * block_size
        X[:, cs] = rng.standard_normal(n)
        for j in range(cs + 1, cs + block_size):
            X[:, j] = rho * X[:, j - 1] + innov_sd * rng.standard_normal(n)

    beta = np.zeros(p)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), n_blocks=n_blocks, rho=rho, snr=snr)


def dgp_ar1_block_white_snr(n, p, support, amplitude=3.0, p_ar=25, p_white=975,
                            n_blocks=5, rho=0.7, snr=2.0, seed=None):
    """M AR(1) blocks (first p_ar columns) plus p_white i.i.d. N(0,1) white-noise
    columns; total width p_ar + p_white. Actives lie in the AR part. SNR-calibrated
    noise. (`p` is the total width supplied by the generic runner; p_ar/p_white
    carry the split.)
    """
    assert p_ar % n_blocks == 0 and p_white >= 1 and abs(rho) < 1 and snr > 0
    rng = np.random.default_rng(seed)
    block_size = p_ar // n_blocks
    innov_sd = np.sqrt(1.0 - rho * rho)

    X_ar = np.zeros((n, p_ar))
    for m in range(n_blocks):
        cs = m * block_size
        X_ar[:, cs] = rng.standard_normal(n)
        for j in range(cs + 1, cs + block_size):
            X_ar[:, j] = rho * X_ar[:, j - 1] + innov_sd * rng.standard_normal(n)
    X_white = rng.standard_normal((n, p_white))
    X = np.hstack([X_ar, X_white])

    beta = np.zeros(p_ar + p_white)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p_ar + p_white, p_ar=p_ar, p_white=p_white, s=len(support),
                n_blocks=n_blocks, rho=rho, snr=snr)


def dgp_ht_block_snr(n, p, support, amplitude=1.0, n_blocks=5, rho=0.8, snr=2.0,
                     nu=3.0, heavy_noise=False, seed=None):
    """Heavy-tailed block DGP: Gaussian AR(1) blocks Z, converted to multivariate
    Student-t(nu) via a Gaussian scale mixture (one chi-squared scale per row,
    shared across columns). Noise is Gaussian, or Student-t(nu) if heavy_noise.
    SNR-calibrated.
    """
    assert p % n_blocks == 0 and abs(rho) < 1 and snr > 0 and nu > 2
    rng = np.random.default_rng(seed)
    block_size = p // n_blocks
    innov_sd = np.sqrt(1.0 - rho * rho)

    Z = np.zeros((n, p))
    for m in range(n_blocks):
        cs = m * block_size
        Z[:, cs] = rng.standard_normal(n)
        for j in range(cs + 1, cs + block_size):
            Z[:, j] = rho * Z[:, j - 1] + innov_sd * rng.standard_normal(n)

    scale = np.sqrt(rng.chisquare(nu, size=n) / nu)   # one scale per observation
    X = Z / scale[:, None]

    beta = np.zeros(p)
    beta[support] = amplitude
    signal = X @ beta
    noise_sigma = np.sqrt(np.var(signal, ddof=1) / snr)
    if heavy_noise:
        eps = rng.standard_t(nu, size=n) * (noise_sigma / np.sqrt(nu / (nu - 2)))
    else:
        eps = rng.standard_normal(n) * noise_sigma
    y = signal + eps

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), n_blocks=n_blocks, rho=rho, snr=snr,
                nu=nu, heavy_noise=heavy_noise)


def dgp_ht_block_white_snr(n, p, support, amplitude=1.0, p_ar=25, p_white=475,
                           n_blocks=5, rho=0.8, snr=2.0, nu=3.0, heavy_noise=False,
                           seed=None):
    """Heavy-tailed AR(1) blocks (p_ar columns) plus a heavy-tailed white block
    (p_white columns), each block scale-mixed to Student-t(nu) with its own
    per-observation chi-squared scale. Actives lie in the AR part. SNR-calibrated.
    """
    assert p_ar % n_blocks == 0 and abs(rho) < 1 and snr > 0 and nu > 2
    rng = np.random.default_rng(seed)
    ptot = p_ar + p_white
    block_size = p_ar // n_blocks
    innov_sd = np.sqrt(1.0 - rho * rho)

    X = np.zeros((n, ptot))
    for m in range(n_blocks):
        cs = m * block_size
        Zb = np.zeros((n, block_size))
        Zb[:, 0] = rng.standard_normal(n)
        for j in range(1, block_size):
            Zb[:, j] = rho * Zb[:, j - 1] + innov_sd * rng.standard_normal(n)
        scale = np.sqrt(rng.chisquare(nu, size=n) / nu)
        X[:, cs:cs + block_size] = Zb / scale[:, None]
    if p_white > 0:
        Zw = rng.standard_normal((n, p_white))
        scale_w = np.sqrt(rng.chisquare(nu, size=n) / nu)
        X[:, p_ar:] = Zw / scale_w[:, None]

    beta = np.zeros(ptot)
    beta[support] = amplitude
    signal = X @ beta
    noise_sigma = np.sqrt(np.var(signal, ddof=1) / snr)
    if heavy_noise:
        eps = rng.standard_t(nu, size=n) * (noise_sigma / np.sqrt(nu / (nu - 2)))
    else:
        eps = rng.standard_normal(n) * noise_sigma
    y = signal + eps

    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=ptot, p_ar=p_ar, p_white=p_white, s=len(support),
                n_blocks=n_blocks, rho=rho, snr=snr, nu=nu, heavy_noise=heavy_noise)


def dgp_groups_toeplitz_leaf(n, p, support, amplitude, groups, rho_levels, phi,
                             snr=2.0, seed=None):
    """Multi-level latent-factor DGP with Toeplitz leaf blocks. `groups` is a list
    of L label-vectors (each length p, fine->coarse); rho_levels is strictly
    decreasing in (0,1). Coarse levels x=1..L-1 add a shared factor per group
    with loading sqrt(rho_ext[x]-rho_ext[x+1]); the finest level adds a Toeplitz(phi)
    block within each leaf; plus idiosyncratic noise. SNR-calibrated.
    """
    assert len(groups) == len(rho_levels) >= 2
    assert all(rho_levels[i] > rho_levels[i + 1] for i in range(len(rho_levels) - 1))
    assert rho_levels[0] < 1.0 and rho_levels[-1] > 0.0 and 0.0 <= phi < 1.0

    rng = np.random.default_rng(seed)
    groups = [np.asarray(g) for g in groups]
    rho_ext = list(rho_levels) + [0.0]
    X = np.zeros((n, p))

    # Coarser latent-factor levels x = 1 .. L-1
    for x in range(1, len(groups)):
        loading = np.sqrt(rho_ext[x] - rho_ext[x + 1])
        for gid in np.unique(groups[x]):
            f = rng.standard_normal(n)
            X[:, groups[x] == gid] += loading * f[:, None]

    # Finest level: Toeplitz(phi) within each leaf block
    leaf_scale = np.sqrt(rho_levels[0] - rho_levels[1])
    for gid in np.unique(groups[0]):
        cols = np.flatnonzero(groups[0] == gid)
        g = len(cols)
        idx = np.arange(g)
        Sigma = phi ** np.abs(idx[:, None] - idx[None, :])
        L = np.linalg.cholesky(Sigma)
        X[:, cols] += leaf_scale * (rng.standard_normal((n, g)) @ L.T)

    # Idiosyncratic noise
    X += np.sqrt(1.0 - rho_levels[0]) * rng.standard_normal((n, p))

    beta = np.zeros(p)
    beta[support] = amplitude
    y = _snr_response(X, beta, snr, rng)
    return dict(X=X, y=y, beta=beta, true_support=np.asarray(support, dtype=np.intp),
                n=n, p=p, s=len(support), snr=snr)
