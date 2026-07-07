# ==============================================================================
# trex_helpers.py
# ==============================================================================
#
# Shared helper functions for TRex single-run demos.
# Mirrors the single-run printing/saving helpers in
# R/trex_selector_methods/trex/trex_sim_utils.R.
#
# ==============================================================================

import numpy as np


def print_single_run_result(scenario_name, dat, sel, tFDR):
    """
    Print a formatted single-run result block.

    Mirrors R .print_single_run().

    Parameters
    ----------
    scenario_name : str
        Label printed in the header (e.g. "Part 1 — i.i.d. Gaussian [TLARS]").
    dat           : dict
        Output of dgp_gauss_snr() with keys X, y, beta, true_support, n, p, s, snr.
    sel           : trex_selector_neo.TRexSelector
        Selector instance after .select() has been called.
    tFDR          : float
        Target FDR level.
    """
    from trex_sim_common import compute_fdp, compute_tpp  # local import avoids circular dep

    true_support = dat["true_support"]          # 0-based np.ndarray
    s            = dat["s"]
    n            = dat["n"]
    p            = dat["p"]

    selected  = list(sel.selected_indices)      # 0-based list
    n_sel     = len(selected)
    true_set  = set(true_support.tolist())
    sel_set   = set(selected)
    n_tp      = len(sel_set & true_set)
    n_fp      = n_sel - n_tp

    fdp_val = compute_fdp(selected, true_support.tolist())
    tpp_val = compute_tpp(selected, true_support.tolist())

    sep  = "=" * 70
    dash = "-" * 70

    print(sep)
    print(f"  {scenario_name}")
    print(dash)
    print(f"  Data:  n = {n},  p = {p},  s = {s},  tFDR = {tFDR:.2f}")
    print(f"         True support (0-based):   {{{', '.join(str(x) for x in true_support)}}}")
    print(f"         Active coefficients:       "
          f"{{{', '.join(f'{x:.3f}' for x in dat['beta'][true_support])}}}")
    print(dash)
    print(f"  Calibration:  T_stop = {sel.T_stop},  L = {sel.L}")
    print(f"  Selection:    {n_sel} selected  |  TP = {n_tp}  FP = {n_fp}")
    print(f"  Rates:        TPP = {tpp_val:.4f}  |  FDP = {fdp_val:.4f}"
          f"  (target tFDR <= {tFDR:.2f})")
    print(dash)

    print("\nAdjusted Relative Occurrences (phi_prime):")
    phi_prime = np.asarray(sel.phi_prime)
    print("  " + "  ".join(f"{x:.6f}" for x in phi_prime))

    print("\nPhi Matrix (phi_mat):")
    print(np.round(np.asarray(sel.phi_mat), 4))

    print("\nEstimated FDP Matrix (fdp_hat_mat):")
    print(np.round(np.asarray(sel.fdp_hat_mat), 4))

    print("\nR Matrix (r_mat):")
    print(np.round(np.asarray(sel.r_mat), 4))

    print("\nVoting Grid (voting_grid):")
    print("  " + "  ".join(f"{x:.4f}" for x in np.asarray(sel.voting_grid)))

    print(sep + "\n")


def save_selection_csv(out_dir, stem, p, phi_prime, selected_0based, true_support_0based):
    """
    Save a per-variable selection CSV.

    Mirrors R .print_and_save_d02() CSV section and C++ save_selection_csv().

    One row per predictor with columns:
        variable_index  (0-based int, matches C++ output)
        phi_prime       (float, 6 d.p.)
        selected        (0 or 1)
        true_positive   (0 or 1)

    Parameters
    ----------
    out_dir              : str        Output directory (created if absent).
    stem                 : str        Base file name without extension.
    p                    : int        Total number of predictors.
    phi_prime            : array-like Shape (p,) adjusted relative occurrences.
    selected_0based      : iterable   0-based selected predictor indices.
    true_support_0based  : iterable   0-based true active predictor indices.
    """
    import os
    import pandas as pd

    os.makedirs(out_dir, exist_ok=True)

    phi = np.asarray(phi_prime, dtype=np.float64)
    sel_set  = set(int(x) for x in selected_0based)
    true_set = set(int(x) for x in true_support_0based)

    df = pd.DataFrame({
        "variable_index": np.arange(p, dtype=np.int64),
        "phi_prime":      np.round(phi, 6),
        "selected":       np.array([1 if i in sel_set  else 0 for i in range(p)], dtype=np.int8),
        "true_positive":  np.array([1 if i in true_set else 0 for i in range(p)], dtype=np.int8),
    })

    csv_path = os.path.join(out_dir, f"{stem}.csv")
    df.to_csv(csv_path, index=False)
    print(f"[Info] CSV saved to: {csv_path}\n")
