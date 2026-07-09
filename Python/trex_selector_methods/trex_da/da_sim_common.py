# ==============================================================================
# da_sim_common.py
# ==============================================================================
#
# Shared infrastructure for the DA-TRex demos: support generators, a DA selector
# builder, a generic parallel Monte Carlo runner, and table/matrix printers.
#
# Mirrors cpp/trex_selector_methods/trex_da/trex_da_sim_common.hpp and the R
# support_generators.R / simulation_utils.R helpers. Support indices are 0-based.
# ==============================================================================

import os
import multiprocessing

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.ml_methods.clustering import LinkageMethod
from trex_selector_neo.utils import compute_fdp, compute_tpp

N_CORES = min(6, os.cpu_count() or 1)

_DA_METHOD = {"AR1": ts.DAMethod.AR1, "NN": ts.DAMethod.NN, "EQUI": ts.DAMethod.EQUI,
              "BT": ts.DAMethod.BT, "PRIOR_GROUPS": ts.DAMethod.PRIOR_GROUPS}
_SOLVER = {"TLARS": ts.SolverTypeForTRex.TLARS,
           "TAFS": ts.SolverTypeForTRex.TAFS,
           "TOMP": ts.SolverTypeForTRex.TOMP}
_LINKAGE = {"single": LinkageMethod.Single, "complete": LinkageMethod.Complete,
            "average": LinkageMethod.Average}


# ==============================================================================
# Support generators (0-based)
# ==============================================================================

def make_support_capped_spread(s, p, max_gap=10):
    """Evenly spaced support with a capped gap: {0, gap, 2*gap, ...}."""
    gap = min(p // s, max_gap)
    return np.arange(s) * gap


def make_support_random(s, p, seed):
    """Uniformly random support without replacement, sorted ascending."""
    rng = np.random.default_rng(seed + 500000)
    return np.sort(rng.choice(p, size=s, replace=False))


def make_support_bt_one_per_block(p, n_blocks):
    """One active predictor per block: the first index of each block."""
    block_size = p // n_blocks
    return np.arange(n_blocks) * block_size


# ==============================================================================
# DA selector
# ==============================================================================

def build_da_selector(method, X, y, tFDR, K, solver="TLARS",
                      rho_thr_DA=None, cor_coef=None, hc_linkage=None,
                      prior_groups=None, rho_grid_labels=None):
    da = ts.TRexDAControlParameter()
    da.method = _DA_METHOD[method]
    if rho_thr_DA is not None:
        da.rho_thr_DA = rho_thr_DA
    if cor_coef is not None:
        da.cor_coef = cor_coef
    if hc_linkage is not None:
        da.hc_linkage = _LINKAGE[hc_linkage]
    if prior_groups is not None:
        da.prior_groups = prior_groups
    if rho_grid_labels is not None:
        da.rho_grid_labels = rho_grid_labels
    tc = ts.TRexControlParameter()
    tc.K = K
    tc.solver_type = _SOLVER[solver]
    sel = ts.TRexDASelector(np.asfortranarray(X), y, tFDR=tFDR,
                            da_control=da, trex_control=tc, seed=-1, verbose=False)
    sel.select()
    return sel


# ==============================================================================
# Monte Carlo runner (parallel; data seed = base + trial, selector seed = -1)
# ==============================================================================

def _mc_worker(args):
    import importlib
    dg = importlib.import_module("dgp_generators")
    (trial_seed, n, p, support, s, amplitude, snr, tFDR, K,
     dgp_name, dgp_extra, da_method, da_extra, solver) = args
    supp = support if support is not None else make_support_random(s, p, trial_seed)
    dat = getattr(dg, dgp_name)(n, p, supp, amplitude=amplitude, snr=snr,
                                seed=trial_seed, **dgp_extra)
    sel = build_da_selector(da_method, dat["X"], dat["y"], tFDR, K, solver, **da_extra)
    selected = list(sel.selected_indices)
    true_support = np.flatnonzero(dat["beta"]).tolist()
    return compute_fdp(selected, true_support), compute_tpp(selected, true_support)


def run_mc(n, p, support, amplitude, snr, tFDR, K, num_MC, seed, dgp_name,
           dgp_extra, da_method, da_extra, s=None, solver="TLARS", n_cores=N_CORES):
    """Run num_MC trials; support=None redraws a random support each trial.
    Returns dict(mean_FDP, mean_TPP, sd_FDP, sd_TPP)."""
    supp = None if support is None else list(np.asarray(support).tolist())
    tasks = [(seed + i, n, p, supp, s, amplitude, snr, tFDR, K,
              dgp_name, dgp_extra, da_method, da_extra, solver) for i in range(num_MC)]
    if n_cores > 1:
        with multiprocessing.get_context("spawn").Pool(n_cores) as pool:
            res = pool.map(_mc_worker, tasks)
    else:
        res = [_mc_worker(t) for t in tasks]

    fdp = np.array([r[0] for r in res])
    tpp = np.array([r[1] for r in res])
    return dict(mean_FDP=fdp.mean(), mean_TPP=tpp.mean(),
                sd_FDP=fdp.std(ddof=1), sd_TPP=tpp.std(ddof=1))


# ==============================================================================
# Printers (mirror R .print_table / .print_matrix)
# ==============================================================================

def print_table(rows, param_col, value_fmt="{:<10.4f}", width=10):
    print(f" {param_col:<{width}} {'mean_FDP':<10} {'mean_TPP':<10} {'sd_FDP':<10} {'sd_TPP':<10}")
    print("-" * 58)
    for r in rows:
        print(" " + value_fmt.format(r[param_col]) + " "
              + f"{r['mean_FDP']:<10.4f} {r['mean_TPP']:<10.4f}"
              + f" {r['sd_FDP']:<10.4f} {r['sd_TPP']:<10.4f}")
    print()


def print_table_multi(rows, key_cols, key_fmts, key_widths):
    """Table with multiple key columns (e.g. Q, p) plus the four metrics.
    key_fmts / key_widths are dicts keyed by column name."""
    metric_cols = ("mean_FDP", "mean_TPP", "sd_FDP", "sd_TPP")
    header = " "
    for k in key_cols:
        header += f"{k:<{key_widths[k]}} "
    for m in metric_cols:
        header += f"{m:<10} "
    print(header.rstrip())
    print("-" * (len(header)))
    for r in rows:
        line = " "
        for k in key_cols:
            line += key_fmts[k].format(r[k]) + " "
        for m in metric_cols:
            line += f"{r[m]:<10.4f} "
        print(line.rstrip())
    print()


def print_matrix(mat, label, row_labels, col_labels):
    print(f"\n {label}:")
    header = f" {'':<10}"
    for cl in col_labels:
        header += f" {cl:<8}"
    print(header)
    print("-" * (len(header) + 2))
    for i, rl in enumerate(row_labels):
        row = f" {rl:<10}"
        for j in range(len(col_labels)):
            row += f" {mat[i, j]:<8.3f}"
        print(row)
    print()
