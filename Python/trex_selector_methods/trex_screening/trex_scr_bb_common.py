# ==============================================================================
# trex_scr_bb_common.py
# ==============================================================================
#
# Shared infrastructure for the Biobank Screen-TRex Python demos (04-05).
#
# Mirrors the C++ demo-internal header
#   cpp/trex_selector_methods/trex_screening/demo_trex_scr_bb_common.hpp
# and Section 5 of the R helper trex_scr_sim_utils.R.
#
# The biobank workflow ("Algorithm 1") screens MANY phenotypes against ONE
# shared design matrix X. Per phenotype it routes adaptively:
#   1. Screen-TRex Ordinary; accept if its estimated FDR lands in the window
#      [lower_bound_FDR, upper_bound_FDR];
#   2. else Screen-TRex Bootstrap-CI;
#   3. else fall back to the classical T-Rex selector at target_FDR_trex.
# The fraction of phenotypes routed to each method is the "Usage %".
#
# Support sets and selected indices are 0-BASED (matching the C++ sources).
# ==============================================================================

import os

import numpy as np
import trex_selector_neo as ts
from trex_selector_neo.utils import compute_fdp, compute_tpp

# Canonical method-name keys returned in each result's `method_used` field.
BIOBANK_METHODS = ["Screen-TRex (ordinary)",
                   "Screen-TRex (bootstrap-CI)",
                   "T-Rex (fallback)"]


# ==============================================================================
# Control factory
# ==============================================================================

def make_biobank_ctrl(use_mmap=False, K=20, R_boot=1000, ci_grid_step=0.001,
                      lower_bound_FDR=0.05, upper_bound_FDR=0.15,
                      target_FDR_trex=0.10):
    """Build a BiobankScreenTRexControl (mirrors makeBiobankControl()).

    BiobankScreenTRexControl nests ScreenTRexControlParameter as
    trex_screen_ctrl, which in turn nests its own trex_ctrl.
    """
    ctrl = ts.BiobankScreenTRexControl()
    ctrl.lower_bound_FDR = lower_bound_FDR
    ctrl.upper_bound_FDR = upper_bound_FDR
    ctrl.target_FDR_trex = target_FDR_trex

    ctrl.trex_screen_ctrl.trex_method = ts.ScreenTRexMethod.TREX
    ctrl.trex_screen_ctrl.R_boot = int(R_boot)
    ctrl.trex_screen_ctrl.ci_grid_step = float(ci_grid_step)

    tc = ctrl.trex_screen_ctrl.trex_ctrl
    tc.K = int(K)
    tc.max_dummy_multiplier = 10
    tc.use_max_T_stop = True
    tc.dummy_distribution = ts.DummyDistribution.normal()
    tc.lloop_strategy = ts.LLoopStrategy.STANDARD
    tc.solver_type = ts.SolverTypeForTRex.TLARS
    tc.use_memory_mapping = bool(use_mmap)
    return ctrl


# ==============================================================================
# Single-phenotype result output
# ==============================================================================

def print_biobank_single(res, true_support):
    """Print a single Biobank Screen-TRex phenotype result."""
    true_support = sorted(int(i) for i in true_support)
    sel = sorted(int(i) for i in res.selected_indices)
    print(f"\n  Method used:       {res.method_used}")
    print(f"  Estimated FDR:     {res.estimated_FDR:.4f}")
    print(f"  a_hat  (ordinary):  {res.estimated_FDR_screen_ordinary:.4f}")
    print(f"  a_hat_C (bootstrap): {res.estimated_FDR_screen_bootstrap:.4f}")
    print(f"  Selected ({len(sel)}):  {' '.join(map(str, sel))}")
    fdp = compute_fdp(sel, true_support)
    tpp = compute_tpp(sel, true_support)
    print(f"  FDP: {fdp:.4f}  |  TPP: {tpp:.4f}")
    print(f"  True support:   {' '.join(map(str, true_support))}\n")


# ==============================================================================
# Multi-phenotype table + summary output
# ==============================================================================

def print_biobank_table(results, true_supports):
    """Print a per-phenotype results table (mirrors print_biobank_phenotype_table)."""
    print("\n" + "=" * 80)
    print("=== Results by Phenotype ===")
    print("=" * 80 + "\n")
    print(f"{'Pheno':>8}{'  Method':<28}{'Est.FDR':>12}{'Selected':>12}"
          f"{'FDP':>12}{'TPP':>12}")
    print("-" * 80)
    for i, r in enumerate(results):
        sel = list(r.selected_indices)
        fdp = compute_fdp(sel, list(true_supports[i]))
        tpp = compute_tpp(sel, list(true_supports[i]))
        print(f"{i:>8}{'  ' + r.method_used:<28}{r.estimated_FDR:>12.4f}"
              f"{len(sel):>12}{fdp:>12.4f}{tpp:>12.4f}")
    print("")


def print_biobank_summary(results):
    """Print multi-phenotype screening summary statistics."""
    n = len(results)
    n_ord = sum(1 for r in results if r.method_used == "Screen-TRex (ordinary)")
    n_boot = sum(1 for r in results if r.method_used == "Screen-TRex (bootstrap-CI)")
    n_trex = sum(1 for r in results if r.method_used == "T-Rex (fallback)")
    avg_fdr = float(np.mean([r.estimated_FDR for r in results]))
    avg_sel = float(np.mean([len(r.selected_indices) for r in results]))
    print("=" * 70)
    print("=== Summary ===")
    print("=" * 70)
    print(f"  Total phenotypes:        {n}")
    print(f"  Screen-TRex (ordinary):  {n_ord}")
    print(f"  Screen-TRex (bootstrap): {n_boot}")
    print(f"  T-Rex (fallback):        {n_trex}")
    print(f"  Avg estimated FDR:       {avg_fdr:.4f}")
    print(f"  Avg selected count:      {avg_sel:.2f}\n")


# ==============================================================================
# MC results: save + print (TXT + CSV, with Usage % rows)
# ==============================================================================

def _bb_row(name, metric, vec, dp, name_w=31, metric_w=12, snr_w=10, col_w=10):
    s = f"{name:<{name_w}}{metric:<{metric_w}}{'':<{snr_w}}"
    s += "".join(f"{v:>{col_w}.{dp}f}" for v in vec)
    return s


def save_and_print_biobank_mc(num_MC, out_dir, file_stem, snr_values,
                              method_names, fdp_map, tpp_map, est_fdr_map,
                              usage_map):
    """Save + print a Biobank MC table with per-method Usage (%) rows."""
    os.makedirs(out_dir, exist_ok=True)
    lines = []
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"=== Biobank Screen-TRex Results (averaged over {num_MC} "
                 f"Monte Carlo runs) ===")
    lines.append("=" * 70)
    lines.append("")

    hdr = f"{'Method':<31}{'Metric':<12}{'SNR':<10}"
    hdr += "".join(f"{s:>10.2f}" for s in snr_values)
    lines.append(hdr)
    lines.append("-" * len(hdr))

    for nm in method_names:
        usage_pct = [u * 100.0 for u in usage_map[nm]]
        lines.append(_bb_row(nm, "Usage (%)", usage_pct,        1))
        lines.append(_bb_row("", "FDR",       fdp_map[nm],      4))
        lines.append(_bb_row("", "TPR",       tpp_map[nm],      4))
        lines.append(_bb_row("", "Est. FDR",  est_fdr_map[nm],  4))
        lines.append("")

    print("\n".join(lines))
    txt_path = os.path.join(out_dir, file_stem + ".txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[Info] TXT results saved to: {txt_path}")

    csv_path = os.path.join(out_dir, file_stem + ".csv")
    with open(csv_path, "w") as fh:
        fh.write("method,metric,snr,value\n")
        for nm in method_names:
            for i, snr in enumerate(snr_values):
                fh.write(f"{nm},Usage,{snr},{usage_map[nm][i]:.6f}\n")
                fh.write(f"{nm},FDR,{snr},{fdp_map[nm][i]:.6f}\n")
                fh.write(f"{nm},TPR,{snr},{tpp_map[nm][i]:.6f}\n")
                fh.write(f"{nm},Est. FDR,{snr},{est_fdr_map[nm][i]:.6f}\n")
    print(f"[Info] CSV results saved to: {csv_path}\n")

# ==============================================================================
# End of trex_scr_bb_common.py
# ==============================================================================
