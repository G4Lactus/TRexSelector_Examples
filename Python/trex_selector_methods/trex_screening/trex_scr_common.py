# ==============================================================================
# trex_scr_common.py
# ==============================================================================
#
# Shared Monte Carlo infrastructure for the Screen-TRex Python demo suite.
#
# Mirrors the C++ demo-internal headers
#   cpp/trex_selector_methods/trex_screening/trex_screening_dgps.hpp
#   cpp/trex_selector_methods/trex_screening/trex_screening_simulation_utils.hpp
#
# Provides:
#   Section 1 - data-generating processes (i.i.d., AR(1), equi, block-equi),
#   Section 2 - Screen-TRex control-parameter factories (flat + object),
#   Section 3 - ProcessPoolExecutor worker + parallel MC runner,
#   Section 4 - dual console/TXT table + tidy pandas CSV output.
#
# Support sets and selected indices are 0-BASED throughout (matching the C++
# sources; the R suite is 1-based). The C++ demos draw with std::mt19937 and
# this suite uses NumPy's Generator, so seeds are mirrored as *labels* only:
# results are statistically comparable to the C++/R reference, not bit-identical.
#
# Parallelism follows the suite convention (see Python/README.md): trials run in
# a concurrent.futures.ProcessPoolExecutor. Workers are module-level functions
# and receive only picklable flat dicts; the trex_selector_neo enums/controls
# are rebuilt inside each worker, never sent across a process boundary. On macOS
# the spawn start method copies sys.path, so the sys.path insert below lets
# spawned workers import this module and its siblings.
# ==============================================================================

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import trex_selector_neo as tsn
from trex_selector_neo.utils import compute_fdp, compute_tpp

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)


# ==============================================================================
# Section 1 - Data-generating processes (0-based supports)
# ==============================================================================

def make_beta(p, p1, beta_val):
    """Evenly-spaced sparse beta vector (mirrors make_beta() in the cpp DGPs).

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
    """i.i.d. Gaussian DGP (mirrors make_iid_dgp / datagen::SyntheticData).

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
    """AR(1) column-correlated DGP (mirrors make_ar1_dgp)."""
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
    """Equicorrelated DGP (mirrors make_equi_dgp).

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
    """Block-equicorrelated DGP (mirrors make_block_equi_dgp).

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
#
# Solver names and screening methods are resolved to trex_selector_neo enums by
# name via getattr, so the flat arg dicts stay picklable across the process
# boundary and the enums are only ever built inside a worker.
#


def _build_trex_ctrl(d):
    """Build a tsn.TRexControlParameter from a flat, picklable dict.

    Mirrors make_trex_control() / make_trex_control_mmap() / make_trex_ctrl_for()
    in the cpp simulation utils. Rebuilt inside worker processes.
    """
    ctrl = tsn.TRexControlParameter()
    ctrl.K = int(d.get("K", 20))
    ctrl.lloop_strategy = tsn.LLoopStrategy.STANDARD
    ctrl.dummy_distribution = tsn.DummyDistribution.normal()
    ctrl.solver_type = getattr(tsn.SolverTypeForTRex, d.get("solver", "TLARS"))
    ctrl.use_memory_mapping = bool(d.get("use_memory_mapping", False))
    ctrl.tloop_stagnation_stop = bool(d.get("tloop_stagnation_stop", False))
    ctrl.tloop_max_stagnant_steps = int(d.get("tloop_max_stagnant_steps", 10))
    if d.get("rho_afs") is not None:
        ctrl.solver_params.rho_afs = float(d["rho_afs"])
    return ctrl


def _build_screen_ctrl(d):
    """Build a tsn.ScreenTRexControlParameter from a flat dict.

    Mirrors make_screen_control(): trex_method / use_bootstrap_CI / R_boot /
    ci_grid_step / rho_thr_DA / n_blocks. Rebuilt inside worker processes.
    """
    sc = tsn.ScreenTRexControlParameter()
    sc.trex_method = getattr(tsn.ScreenTRexMethod, d.get("trex_method", "TREX"))
    sc.use_bootstrap_CI = bool(d.get("bootstrap", False))
    sc.R_boot = int(d.get("R_boot", 1000))
    sc.ci_grid_step = float(d.get("ci_grid_step", 0.001))
    sc.rho_thr_DA = float(d.get("rho_thr_DA", 0.02))
    sc.n_blocks = int(d.get("screen_n_blocks", 5))
    return sc


def default_scr_methods():
    """Ordinary-vs-Bootstrap method descriptors (mirrors default_methods()).

    Names match the cpp display strings.
    """
    return [
        {"name": "Screen-TRex Ordinary: TLARS",     "trex_method": "TREX", "bootstrap": False},
        {"name": "Screen-TRex Bootstrap-CI: TLARS", "trex_method": "TREX", "bootstrap": True},
    ]


# ==============================================================================
# Section 3 - ProcessPoolExecutor worker + parallel MC runner
# ==============================================================================

def _build_screen_data(args):
    """Generate one Monte Carlo dataset from a flat args dict (inside a worker).

    Supports the four DGP kinds used across the suite. For the i.i.d. DGP the
    active support and (optionally) coefficients are drawn from RNG streams
    offset from the trial seed (by 500000 / 600000), matching the cpp make_iid_dgp
    isolation of support/coefficient draws from the design/noise draws.
    """
    kind = args["dgp_kind"]
    n, p = args["n"], args["p"]
    snr = args["snr"]
    seed = args["trial_seed"]

    if kind == "iid":
        support_size = args["support_size"]
        support = args.get("fixed_support")
        if support is None:
            rng_sup = np.random.default_rng(seed + 500000)          # isolated support RNG
            support = rng_sup.choice(p, size=support_size, replace=False)
        else:
            support = np.asarray(support)
        if args.get("rnd_coef", False):
            rng_coef = np.random.default_rng(seed + 600000)         # isolated coefficient RNG
            coefs = rng_coef.standard_normal(len(support))
        else:
            coefs = np.full(len(support), 1.0)
        return dgp_iid_snr(n, p, support, coefs, snr, seed)

    p1, rho, beta_val = args["p1"], args["rho"], args["beta_val"]
    if kind == "ar1":
        return dgp_ar1(n, p, p1, rho, snr, beta_val, seed)
    if kind == "equi":
        return dgp_equi(n, p, p1, rho, snr, beta_val, seed)
    if kind == "block_equi":
        return dgp_block_equi(n, p, p1, rho, args["dgp_n_blocks"], snr, beta_val, seed)
    raise ValueError(f"unknown dgp_kind: {kind!r}")


def _screen_trial_worker(args):
    """Top-level ProcessPoolExecutor worker: one Screen-TRex Monte Carlo trial.

    Must remain module-level so pickle can resolve it as
    ('trex_scr_common', '_screen_trial_worker'). Builds its own data and controls,
    runs the selector with seed = -1 (hardware entropy per trial, as in the cpp
    MC loop), and returns picklable per-trial metrics.
    """
    import trex_selector_neo as _tsn
    from trex_selector_neo.utils import compute_fdp as _fdp, compute_tpp as _tpp

    dat = _build_screen_data(args)
    trex_ctrl = _build_trex_ctrl(args)
    screen_ctrl = _build_screen_ctrl(args)

    sel = _tsn.TRexScreeningSelector(
        dat["X"], dat["y"], screen_control=screen_ctrl,
        trex_control=trex_ctrl, seed=-1, verbose=False)
    res = sel.select()

    selected = list(sel.selected_indices)
    tsup = list(dat["true_support"])
    return {"fdp": _fdp(selected, tsup),
            "tpp": _tpp(selected, tsup),
            "est": float(res.estimated_FDR)}


def run_mc_screen(num_MC, label, dgp_args, trex_args, screen_args,
                  base_seed, max_workers):
    """Run num_MC Screen-TRex trials for one grid point in parallel.

    Mirrors run_mc_trials_screen_trex(): trial mc uses seed base_seed + mc; each
    trial runs an independent selector with seed = -1. Returns
    dict(avg_fdr, avg_tpr, avg_est_fdr, sd_fdr, sd_tpr).

    dgp_args / trex_args / screen_args are flat, picklable dicts merged (with the
    per-trial seed) into the worker arg dict.
    """
    args_list = [
        {**dgp_args, **trex_args, **screen_args, "trial_seed": base_seed + mc}
        for mc in range(num_MC)
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_screen_trial_worker, args_list))

    fdp = np.array([r["fdp"] for r in results])
    tpp = np.array([r["tpp"] for r in results])
    est = np.array([r["est"] for r in results])
    out = {
        "avg_fdr": float(fdp.mean()),
        "avg_tpr": float(tpp.mean()),
        "avg_est_fdr": float(est.mean()),
        "sd_fdr": float(fdp.std(ddof=1)) if num_MC > 1 else 0.0,
        "sd_tpr": float(tpp.std(ddof=1)) if num_MC > 1 else 0.0,
    }
    print(f"  {label} - done.  TPR={out['avg_tpr']:.4f}  "
          f"FDR={out['avg_fdr']:.4f}  Est.FDR={out['avg_est_fdr']:.4f}")
    return out


# ==============================================================================
# Section 4 - Dual console/TXT table + tidy pandas CSV output
# ==============================================================================
#
# Layout mirrors save_and_print_mc_results() in trex_screening_simulation_utils.hpp:
# the method name gets its own line, metric rows are indented beneath it, sweep
# values are right-aligned so the columns stay put regardless of name length.
#
_INDENT_W = 2
_METRIC_W = 23
_COL_W = 10


def save_and_print_scr_mc(num_MC, out_dir, file_stem, x_values, method_names,
                          fdr_map, tpr_map, est_fdr_map, sweep_label="SNR"):
    """Save + print a sweep-variable MC table (mirrors save_and_print_mc_results).

    Writes a human-readable TXT table (also echoed to console) and a tidy
    long-format CSV (method, metric, <sweep_label>, value) into out_dir via pandas.
    """
    os.makedirs(out_dir, exist_ok=True)
    lines = []
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"=== Screen-TRex Results (averaged over {num_MC} Monte Carlo runs) ===")
    lines.append("=" * 70)
    lines.append("")

    hdr = " " * _INDENT_W + f"{sweep_label:<{_METRIC_W}}"
    hdr += "".join(f"{x:>{_COL_W}.2f}" for x in x_values)
    lines.append(hdr)
    sep_w = _INDENT_W + _METRIC_W + _COL_W * len(x_values)
    lines.append("-" * sep_w)

    def _metric_row(label, vec):
        s = " " * _INDENT_W + f"{label:<{_METRIC_W}}"
        s += "".join(f"{v:>{_COL_W}.4f}" for v in vec)
        return s

    for nm in method_names:
        lines.append("")
        lines.append(nm)                       # method name on its own line
        lines.append(_metric_row("FDR",      fdr_map[nm]))
        lines.append(_metric_row("TPR",      tpr_map[nm]))
        lines.append(_metric_row("Est. FDR", est_fdr_map[nm]))
    lines.append("")

    print("\n".join(lines))
    txt_path = os.path.join(out_dir, file_stem + ".txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[Info] TXT results saved to: {txt_path}")

    rows = []
    for nm in method_names:
        for i, x in enumerate(x_values):
            rows.append({"method": nm, "metric": "FDR",      sweep_label: x, "value": fdr_map[nm][i]})
            rows.append({"method": nm, "metric": "TPR",      sweep_label: x, "value": tpr_map[nm][i]})
            rows.append({"method": nm, "metric": "Est. FDR", sweep_label: x, "value": est_fdr_map[nm][i]})
    csv_path = os.path.join(out_dir, file_stem + ".csv")
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"[Info] CSV results saved to: {csv_path}\n")

# ==============================================================================
# End of trex_scr_common.py
# ==============================================================================
