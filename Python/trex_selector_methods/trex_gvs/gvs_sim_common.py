# ==============================================================================
# gvs_sim_common.py
# ==============================================================================
#
# Shared infrastructure for the T-Rex+GVS demos: a GVS selector builder, a
# generic parallel Monte Carlo runner over the trex_gvs_dgps generators, and
# table / matrix / single-run printers.
#
# Mirrors R/trex_selector_methods/simulation_utils.R (the .gvs_select adapter,
# the .run_mc_* GVS wrappers, and .print_table / .print_matrix) and the C++
# trex_gvs_simulation_utils.hpp. Support / group indices are 0-based.
# ==============================================================================

import os
import multiprocessing

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.ml_methods.clustering import LinkageMethod
from trex_selector_neo.utils import compute_fdp, compute_tpp

N_CORES = min(6, os.cpu_count() or 1)

_GVS_TYPE = {"EN": ts.GVSType.EN, "IEN": ts.GVSType.IEN}
_EN_SOLVER = {"TENET": ts.ENSolverType.TENET, "TENET_AUG": ts.ENSolverType.TENET_AUG}
_LAMBDA2 = {"CV_1SE_CCD": ts.LambdaSelectionMethod.CV_1SE_CCD,
            "CV_1SE_SVD": ts.LambdaSelectionMethod.CV_1SE_SVD,
            "CV_MIN_CCD": ts.LambdaSelectionMethod.CV_MIN_CCD,
            "CV_MIN_SVD": ts.LambdaSelectionMethod.CV_MIN_SVD}
_LINKAGE = {"single": LinkageMethod.Single, "complete": LinkageMethod.Complete,
            "average": LinkageMethod.Average, "wpgma": LinkageMethod.WPGMA,
            "ward": LinkageMethod.Ward, "centroid": LinkageMethod.Centroid,
            "median": LinkageMethod.Median}

_HC_LABEL = {"single": "Single", "complete": "Complete", "average": "Average",
             "wpgma": "WPGMA", "ward": "Ward", "centroid": "Centroid",
             "median": "Median"}


def cap_hc(hc_dist):
    """Capitalised display name for a linkage key (e.g. 'single' -> 'Single')."""
    return _HC_LABEL.get(hc_dist, hc_dist)


# ==============================================================================
# GVS selector
# ==============================================================================

def build_gvs_selector(X, y, tFDR, K, gvs_type, corr_max, hc_linkage="single",
                       en_solver="TENET", lambda2_method="CV_1SE_CCD",
                       groups=None, cv_seed=None, return_result=False):
    """Build and run a TRexGVSSelector with the given group-selection knobs.

    Centralises the selector construction shared by every MC runner and demo,
    so solver / linkage / lambda_2 options are wired consistently. `groups`
    supplies oracle group labels (a length-p int vector of contiguous 0-based
    cluster IDs, passed as prior_groups); None (default) discovers groups from
    the data via HAC.

    Returns the selector after select() has been called. With
    `return_result=True`, returns the tuple ``(sel, res)`` where ``res`` is the
    ``GVSSelectionResult`` returned by ``select()`` — it exposes the grouping
    diagnostics ``max_clusters`` (number of clusters formed) and ``groups_vec``
    (per-variable cluster labels) that the selector object itself does not carry.
    """
    gc = ts.TRexGVSControlParameter()
    gc.gvs_type = _GVS_TYPE[gvs_type]
    gc.corr_max = corr_max
    gc.hc_linkage = _LINKAGE[hc_linkage]
    gc.en_solver = _EN_SOLVER[en_solver]
    gc.lambda2_method = _LAMBDA2[lambda2_method]
    if cv_seed is not None:
        gc.cv_seed = cv_seed
    if groups is not None:
        gc.prior_groups = list(np.asarray(groups).astype(int).tolist())

    tc = ts.TRexControlParameter()
    tc.K = K
    tc.solver_type = ts.SolverTypeForTRex.TLARS
    gc.trex_ctrl = tc

    sel = ts.TRexGVSSelector(np.asfortranarray(X), y, tFDR=tFDR,
                             gvs_control=gc, seed=-1, verbose=False)
    res = sel.select()
    return (sel, res) if return_result else sel


# ==============================================================================
# Monte Carlo runner (parallel; data seed = base + trial, selector seed = -1)
# ==============================================================================

def _mc_worker(args):
    import importlib
    dg = importlib.import_module("trex_gvs_dgps")
    (trial_seed, dgp_name, dgp_kwargs, tFDR, K, gvs_type, corr_max,
     hc_linkage, en_solver, lambda2_method) = args
    dat = getattr(dg, dgp_name)(seed=trial_seed, **dgp_kwargs)
    sel = build_gvs_selector(dat["X"], dat["y"], tFDR, K, gvs_type, corr_max,
                             hc_linkage, en_solver, lambda2_method)
    selected = list(sel.selected_indices)
    true_support = np.flatnonzero(dat["beta"]).tolist()
    return compute_fdp(selected, true_support), compute_tpp(selected, true_support)


def run_mc_gvs(dgp_name, dgp_kwargs, tFDR, K, num_MC, seed, gvs_type, corr_max,
               hc_linkage="single", en_solver="TENET",
               lambda2_method="CV_1SE_CCD", n_cores=N_CORES):
    """Run num_MC trials of `dgp_name` (a trex_gvs_dgps generator) + trex+GVS.

    Trial k uses data seed = seed + k; the selector seed is fixed at -1. All
    DGP-specific parameters (n, p, snr, sd_x / rho / ar_coef, ...) are passed via
    `dgp_kwargs`. Returns dict(mean_FDP, mean_TPP, sd_FDP, sd_TPP).
    """
    tasks = [(seed + i, dgp_name, dgp_kwargs, tFDR, K, gvs_type, corr_max,
              hc_linkage, en_solver, lambda2_method) for i in range(num_MC)]
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
# Printers (mirror R .print_table / .print_matrix / .print_result)
# ==============================================================================

def print_table(rows, param_col, value_fmt="{:<10.4f}", width=10):
    print(f" {param_col:<{width}} {'mean_FDP':<10} {'mean_TPP':<10}"
          f" {'sd_FDP':<10} {'sd_TPP':<10}")
    print("-" * 58)
    for r in rows:
        print(" " + value_fmt.format(r[param_col]) + " "
              + f"{r['mean_FDP']:<10.4f} {r['mean_TPP']:<10.4f}"
              + f" {r['sd_FDP']:<10.4f} {r['sd_TPP']:<10.4f}")
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


def print_single_run_result(scenario_name, dat, sel, tFDR):
    """Formatted single-run result block, mirroring the R demos' .print_result.
    Prints the sd_x line only when the DGP exposes it (demos 01-04)."""
    ts_idx = np.flatnonzero(dat["beta"]).tolist()
    sel_set = list(sel.selected_indices)
    n_sel = len(sel_set)
    n_tp = len(set(sel_set) & set(ts_idx))
    n_fp = n_sel - n_tp

    print("=" * 70)
    print(f"  {scenario_name}")
    print("-" * 70)
    print(f" Data: n = {dat['n']}, p = {dat['p']}, tFDR = {tFDR:.2f}, s = {dat['s']}")
    if "sd_x" in dat:
        print(f"       SNR = {dat['snr']:.2f},  sigma_y = {dat['sigma_y']:.4f},"
              f"  sd_x = {dat['sd_x']:.4f}")
    else:
        print(f"       SNR = {dat['snr']:.2f},  sigma_y = {dat['sigma_y']:.4f}")
    print("-" * 70)
    print(f"  Calibration:  T_stop = {sel.T_stop},  L = {sel.L}")
    print(f"  Selection:    {n_sel} selected  |  TP = {n_tp}  FP = {n_fp}")
    print(f"  Rates:        TPP = {compute_tpp(sel_set, ts_idx):.3f}  |"
          f"  FDP = {compute_fdp(sel_set, ts_idx):.3f}  (target tFDR <= {tFDR:.2f})")
    print("=" * 70 + "\n")
