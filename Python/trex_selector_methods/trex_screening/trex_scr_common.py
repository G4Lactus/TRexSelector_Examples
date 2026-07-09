# ==============================================================================
# trex_scr_common.py
# ==============================================================================
#
# Shared infrastructure for the Screen-TRex Python demonstration suite.
#
# Mirrors the C++ demo-internal header
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_common.hpp
# and the Sections 1-4 of the R helper
#   R/trex_selector_methods/trex_screening/trex_scr_sim_utils.R
#
# Provides:
#   Section 1 - data-generating processes (i.i.d., AR(1), equi, block-equi),
#   Section 2 - Screen-TRex control-parameter factories,
#   Section 3 - single-run result printer,
#   Section 4 - Monte Carlo runner + dual TXT/CSV table output.
#
# Support sets and selected indices are 0-BASED throughout (matching the C++
# sources; the R suite is 1-based). The C++ demos draw with std::mt19937 and
# this suite uses NumPy's Generator, so seeds are mirrored as *labels* only:
# results are statistically comparable to the C++/R reference, not bit-identical.
# ==============================================================================

import os

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.utils import compute_fdp, compute_tpp

_SOLVER = {"TLARS": ts.SolverTypeForTRex.TLARS,
           "TAFS":  ts.SolverTypeForTRex.TAFS,
           "TOMP":  ts.SolverTypeForTRex.TOMP}

_SCREEN_METHOD = {"TREX":               ts.ScreenTRexMethod.TREX,
                  "TREX_DA_AR1":        ts.ScreenTRexMethod.TREX_DA_AR1,
                  "TREX_DA_EQUI":       ts.ScreenTRexMethod.TREX_DA_EQUI,
                  "TREX_DA_BLOCK_EQUI": ts.ScreenTRexMethod.TREX_DA_BLOCK_EQUI}


# ==============================================================================
# Section 1 - Data-generating processes (0-based supports)
# ==============================================================================

def make_beta(p, p1, beta_val):
    """Evenly-spaced sparse beta vector (mirrors make_beta() in cpp demo 03).

    Active indices are placed at round(i * p / (p1 + 2)) for i = 1..p1 (0-based).
    """
    beta = np.zeros(p)
    for i in range(1, p1 + 1):
        idx = int(round(i * p / (p1 + 2)))
        beta[idx] = beta_val
    return beta


def _make_response(X, beta, snr, rng):
    """SNR-calibrated response y = X beta + eps, mean-centred.

    sigma = sqrt(||X beta||^2 / (n * snr)); mirrors make_response().
    """
    n = X.shape[0]
    signal = X @ beta
    sigma = np.sqrt(float(signal @ signal) / (n * snr))
    y = signal + rng.standard_normal(n) * sigma
    return y - y.mean()


def _standardise_cols(X):
    """Standardise columns to mean 0 and unit variance (n - 1 denominator)."""
    X = X - X.mean(axis=0, keepdims=True)
    sds = np.sqrt((X ** 2).sum(axis=0) / (X.shape[0] - 1))
    sds[sds < 1e-12] = 1.0
    return X / sds


def dgp_iid_snr(n, p, true_support, coefs, snr, seed):
    """i.i.d. Gaussian DGP (mirrors datagen::SyntheticData in demos 01/02/06).

    X ~ N(0, 1) i.i.d.; y = X beta + eps with SNR-calibrated noise.
    Returns dict(X, y, true_support) with a sorted 0-based support.
    """
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, p))
    beta = np.zeros(p)
    beta[true_support] = coefs
    y = _make_response(X, beta, snr, rng)
    return {"X": X, "y": y, "true_support": np.sort(np.asarray(true_support))}


def dgp_ar1(n, p, p1, rho, snr, beta_val, seed):
    """AR(1) column-correlated DGP (mirrors dgp_ar1 in cpp demo 03)."""
    rng = np.random.default_rng(seed)
    innov_sd = np.sqrt(1 - rho ** 2)
    X = np.zeros((n, p))
    X[:, 0] = rng.standard_normal(n)
    for j in range(1, p):
        X[:, j] = rho * X[:, j - 1] + innov_sd * rng.standard_normal(n)
    X = _standardise_cols(X)
    beta = make_beta(p, p1, beta_val)
    y = _make_response(X, beta, snr, rng)
    return {"X": X, "y": y, "true_support": np.where(beta != 0)[0]}


def dgp_equi(n, p, p1, rho, snr, beta_val, seed):
    """Equicorrelated DGP (mirrors dgp_equi in cpp demo 03).

    x_j = sqrt(rho) * f + sqrt(1 - rho) * w_j, one shared factor f per row.
    """
    rng = np.random.default_rng(seed)
    load_f = np.sqrt(rho)
    load_eta = np.sqrt(1 - rho)
    f = rng.standard_normal(n)
    X = load_f * f[:, None] + load_eta * rng.standard_normal((n, p))
    X = _standardise_cols(X)
    beta = make_beta(p, p1, beta_val)
    y = _make_response(X, beta, snr, rng)
    return {"X": X, "y": y, "true_support": np.where(beta != 0)[0]}


def dgp_block_equi(n, p, p1, rho, n_blocks, snr, beta_val, seed):
    """Block-equicorrelated DGP (mirrors dgp_block_equi in cpp demo 03).

    Contiguous blocks, each with its own shared factor.
    """
    rng = np.random.default_rng(seed)
    load_f = np.sqrt(rho)
    load_eta = np.sqrt(1 - rho)
    base_size = p // n_blocks
    remainder = p % n_blocks
    X = np.zeros((n, p))
    col = 0
    for k in range(n_blocks):
        bk = base_size + (1 if k < remainder else 0)
        z_k = rng.standard_normal(n)                    # shared factor for block k
        X[:, col:col + bk] = (load_f * z_k[:, None]
                              + load_eta * rng.standard_normal((n, bk)))
        col += bk
    X = _standardise_cols(X)
    beta = make_beta(p, p1, beta_val)
    y = _make_response(X, beta, snr, rng)
    return {"X": X, "y": y, "true_support": np.where(beta != 0)[0]}


# ==============================================================================
# Section 2 - Control-parameter factories
# ==============================================================================

def make_scr_trex_ctrl(solver="TLARS", use_memory_mapping=False, **overrides):
    """Default TRexControlParameter for Screen-TRex (mirrors make_scr_trex_ctrl).

    K = 20, max_dummy_multiplier = 10, use_max_T_stop = True, STANDARD L-loop,
    Normal dummy distribution. Extra keyword overrides are applied to the
    control (rho_afs is routed to solver_params; tloop_stagnation_stop and
    tloop_max_stagnant_steps set the stagnation stop used by TAFS/TOMP).
    """
    ctrl = ts.TRexControlParameter()
    ctrl.K = int(overrides.pop("K", 20))
    ctrl.max_dummy_multiplier = int(overrides.pop("max_dummy_multiplier", 10))
    ctrl.use_max_T_stop = bool(overrides.pop("use_max_T_stop", True))
    ctrl.lloop_strategy = ts.LLoopStrategy.STANDARD
    ctrl.dummy_distribution = ts.DummyDistribution.normal()
    ctrl.solver_type = _SOLVER[solver]
    ctrl.use_memory_mapping = bool(use_memory_mapping)
    ctrl.tloop_stagnation_stop = bool(overrides.pop("tloop_stagnation_stop", False))
    ctrl.tloop_max_stagnant_steps = int(overrides.pop("tloop_max_stagnant_steps", 10))
    if "rho_afs" in overrides:
        ctrl.solver_params.rho_afs = float(overrides.pop("rho_afs"))
    if overrides:
        raise TypeError(f"unexpected control overrides: {sorted(overrides)}")
    return ctrl


def make_screen_ctrl(trex_method="TREX", bootstrap=False, R_boot=1000,
                     ci_grid_step=0.001):
    """Default ScreenTRexControlParameter (mirrors make_screen_control())."""
    sc = ts.ScreenTRexControlParameter()
    sc.trex_method = _SCREEN_METHOD[trex_method]
    sc.use_bootstrap_CI = bool(bootstrap)
    sc.R_boot = int(R_boot)
    sc.ci_grid_step = float(ci_grid_step)
    return sc


def default_scr_methods():
    """Ordinary-vs-Bootstrap method descriptors (mirrors default_methods())."""
    return [
        {"name": "Screen-TRex Ordinary",  "trex_method": "TREX", "bootstrap": False},
        {"name": "Screen-TRex Bootstrap", "trex_method": "TREX", "bootstrap": True},
    ]


# ==============================================================================
# Section 3 - Single-run result printer
# ==============================================================================

def print_scr_result(sel, res, true_support):
    """Print a single Screen-TRex selection block (mirrors print_screen_result)."""
    selected = np.sort(sel.selected_indices).tolist()
    true_support = sorted(int(i) for i in true_support)
    print(f"\n  Estimated FDR:  {res.estimated_FDR:.4f}")
    if abs(res.estimated_correlation) > 0:
        print(f"  Estimated rho:  {res.estimated_correlation:.4f}")
    if res.used_bootstrap:
        print("  Bootstrap:      yes")
        print(f"  Conf. level:    {res.confidence_level:.4f}")
    print(f"  Selected ({len(selected)}):  {' '.join(map(str, selected))}")
    fdp = compute_fdp(selected, true_support)
    tpp = compute_tpp(selected, true_support)
    print(f"  FDP: {fdp:.4f}  |  TPP: {tpp:.4f}")
    print(f"  True support:   {' '.join(map(str, true_support))}\n")
    return {"fdp": fdp, "tpp": tpp}


# ==============================================================================
# Section 4 - Monte Carlo runner + table output
# ==============================================================================

def run_mc_screen(num_MC, label, make_data, trex_ctrl, screen_args, base_seed):
    """Run num_MC Screen-TRex trials for one grid point (mirrors run_mc_screen).

    make_data(seed) -> dict(X, y, true_support); the selector uses seed = -1
    (hardware entropy per trial, matching the C++ MC loop). Returns
    dict(avg_fdr, avg_tpr, avg_est_fdr).
    """
    fdp = np.empty(num_MC)
    tpp = np.empty(num_MC)
    est = np.empty(num_MC)
    for mc in range(num_MC):
        dat = make_data(base_seed + mc)
        sc = make_screen_ctrl(trex_method=screen_args["trex_method"],
                              bootstrap=screen_args["bootstrap"],
                              R_boot=screen_args.get("R_boot", 1000))
        sel = ts.TRexScreeningSelector(dat["X"], dat["y"], screen_control=sc,
                                       trex_control=trex_ctrl, seed=-1, verbose=False)
        res = sel.select()
        selected = list(sel.selected_indices)
        tsup = dat["true_support"].tolist()
        fdp[mc] = compute_fdp(selected, tsup)
        tpp[mc] = compute_tpp(selected, tsup)
        est[mc] = res.estimated_FDR
    avg = (float(fdp.mean()), float(tpp.mean()), float(est.mean()))
    print(f"  {label} -- done.  TPR={avg[1]:.4f}  FDR={avg[0]:.4f}  Est.FDR={avg[2]:.4f}")
    return {"avg_fdr": avg[0], "avg_tpr": avg[1], "avg_est_fdr": avg[2]}


def _fmt_row(name, metric, vec, name_w=31, metric_w=10, sweep_w=10, col_w=10):
    s = f"{name:<{name_w}}{metric:<{metric_w}}{'':<{sweep_w}}"
    s += "".join(f"{v:>{col_w}.4f}" for v in vec)
    return s


def save_and_print_scr_mc(num_MC, out_dir, file_stem, x_values, method_names,
                          fdr_map, tpr_map, est_fdr_map, sweep_label="SNR"):
    """Save + print a sweep-variable MC table (mirrors save_and_print_mc_results).

    Writes a human-readable TXT table (also echoed to console) and a tidy
    long-format CSV (method, metric, <sweep_label>, value) into out_dir.
    """
    os.makedirs(out_dir, exist_ok=True)
    lines = []
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"=== Screen-TRex Results (averaged over {num_MC} Monte Carlo runs) ===")
    lines.append("=" * 70)
    lines.append("")

    hdr = f"{'Method':<31}{'Metric':<10}{sweep_label:<10}"
    hdr += "".join(f"{x:>10.2f}" for x in x_values)
    lines.append(hdr)
    lines.append("-" * len(hdr))

    for nm in method_names:
        lines.append(_fmt_row(nm, "FDR",      fdr_map[nm]))
        lines.append(_fmt_row("", "TPR",      tpr_map[nm]))
        lines.append(_fmt_row("", "Est. FDR", est_fdr_map[nm]))
        lines.append("")

    print("\n".join(lines))
    txt_path = os.path.join(out_dir, file_stem + ".txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[Info] TXT results saved to: {txt_path}")

    csv_path = os.path.join(out_dir, file_stem + ".csv")
    with open(csv_path, "w") as fh:
        fh.write(f"method,metric,{sweep_label},value\n")
        for nm in method_names:
            for i, x in enumerate(x_values):
                fh.write(f"{nm},FDR,{x},{fdr_map[nm][i]:.6f}\n")
                fh.write(f"{nm},TPR,{x},{tpr_map[nm][i]:.6f}\n")
                fh.write(f"{nm},Est. FDR,{x},{est_fdr_map[nm][i]:.6f}\n")
    print(f"[Info] CSV results saved to: {csv_path}\n")

# ==============================================================================
# End of trex_scr_common.py
# ==============================================================================
