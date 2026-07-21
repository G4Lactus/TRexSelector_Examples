# ==============================================================================
# demo_trex_gvs_01_mc_sim_hastie_en_blocks.py
# ==============================================================================
#
# MC simulation: T-Rex+GVS on the Hastie DGP.
# This scenario evaluates grouped-signal recovery under strong equicorrelation,
# not sparse within-group support identification. Within-group correlation is
# controlled through sd_x, with rho = 1 / (1 + sd_x^2).
#
# Mirrors cpp/trex_selector_methods/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks.
#
# Three parts, each in its own function:
#   sim_hastie_part_1() — 2-D SNR x rho grid
#   sim_hastie_part_2() — lambda2_method comparison: CV_1SE_SVD / CV_1SE_CCD
#   sim_hastie_part_3() — hc_linkage comparison: Single / Complete / Average / WPGMA
#
# DGP: three equicorrelated groups of 50 variables each, all active (s = 150),
#      plus 350 null white-noise variables; p = 500 and n = 200.
# Selector: K=20, tFDR=0.1, corr_max=0.98, lambda2=CV_1SE_CCD.
# Methods compared: TENET, TENET_AUG, TIENET_AUG.
# MC: 200 trials per grid point (the full run takes hours — reduce num_MC and
#     the grids for a quick smoke test). Optional first CLI argument sets the
#     number of parallel worker processes (default min(6, cpu_count)).
# Indices are 0-based (Python convention); the R counterpart is 1-based.
# ==============================================================================

import math
import os
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from gvs_sim_common import (N_CORES, fmt_num, print_mc_matrix,
                            print_section_header, run_gvs_2d_sweep,
                            run_gvs_mc_trials, save_mc_1d_method_tables,
                            save_mc_2d_tables)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES

OUT_DIR = os.path.join(_THIS_DIR, "simulation_results", "data")

# Simulation parameters (mirrors the cpp GVSSimConfig block in main()).
CFG = dict(
    n=200, p=500,
    sd_x=0.1,           # sqrt(0.01) => rho ~ 0.99 (Parts fix rho locally)
    tFDR=0.1, K=20,
    num_MC=200,
    base_seed=2026,
    corr_max=0.98,
    snr=2.0,            # unused here; Parts 2-3 fix SNR = 2.0 locally
)

SCENARIO_TAG = "Hastie"

RUN_PART_1 = True
RUN_PART_2 = True
RUN_PART_3 = True


# ==============================================================================
#  Part 1 — 2-D SNR x rho sweep
# ==============================================================================

def sim_hastie_part_1(cfg, tag):
    """2-D MC sweep over SNR and rho: TENET / TENET_AUG / TIENET_AUG."""
    snr_grid_2d = [0.2, 0.5, 1.0, 2.0, 5.0]
    rho_grid_2d = [0.30, 0.50, 0.70, 0.90, 0.95, 0.99]

    def cell_kwargs(snr, rho):
        sd_x = math.sqrt((1.0 - rho) / rho)
        return dict(n=cfg["n"], p=cfg["p"], snr=snr, sd_x=sd_x)

    print_section_header(f"Part 1/3: 2-D SNR x rho Sweep  |  {tag}")
    print(f"  n={cfg['n']}  p={cfg['p']}  s=150  MC={cfg['num_MC']}"
          f"  tFDR={fmt_num(cfg['tFDR'])}  K={cfg['K']}"
          f"  corr_max={fmt_num(cfg['corr_max'])}\n"
          "  methods: TENET / TENET_AUG / TIENET_AUG  |  lambda2=CV_1SE_CCD"
          "  linkage=Single\n"
          "  snr_grid: {" + ", ".join(fmt_num(s) for s in snr_grid_2d) + "}\n"
          "  rho_grid: {" + ", ".join(fmt_num(r) for r in rho_grid_2d) + "}\n")

    en_2d = run_gvs_2d_sweep("dgp_hastie_snr", cell_kwargs, snr_grid_2d,
                             rho_grid_2d, cfg, "EN", "TENET", "TENET",
                             max_workers=N_WORKERS)
    en_aug_2d = run_gvs_2d_sweep("dgp_hastie_snr", cell_kwargs, snr_grid_2d,
                                 rho_grid_2d, cfg, "EN", "TENET_AUG",
                                 "TENET_AUG", max_workers=N_WORKERS)
    ien_2d = run_gvs_2d_sweep("dgp_hastie_snr", cell_kwargs, snr_grid_2d,
                              rho_grid_2d, cfg, "IEN", "TENET", "TIENET_AUG",
                              max_workers=N_WORKERS)

    snr_labels = [f"snr={fmt_num(s)}" for s in snr_grid_2d]
    rho_labels = [f"rho={fmt_num(r)}" for r in rho_grid_2d]

    print_mc_matrix("mean_TPP [TENET]", snr_labels, rho_labels, en_2d, True)
    print_mc_matrix("mean_FDP [TENET]", snr_labels, rho_labels, en_2d, False)
    print_mc_matrix("mean_TPP [TENET_AUG]", snr_labels, rho_labels, en_aug_2d, True)
    print_mc_matrix("mean_FDP [TENET_AUG]", snr_labels, rho_labels, en_aug_2d, False)
    print_mc_matrix("mean_TPP [TIENET_AUG]", snr_labels, rho_labels, ien_2d, True)
    print_mc_matrix("mean_FDP [TIENET_AUG]", snr_labels, rho_labels, ien_2d, False)
    save_mc_2d_tables(tag, snr_labels, rho_labels, en_2d, en_aug_2d, ien_2d,
                      cfg, OUT_DIR, f"gvs_{tag}_2d")


# ==============================================================================
#  Part 2 — lambda2_method comparison
# ==============================================================================

def sim_hastie_part_2(cfg, tag):
    """Compare lambda2 selection methods (CV_1SE_SVD vs CV_1SE_CCD) for TENET
    and TIENET_AUG at the fixed operating point SNR = 2.0, rho = 0.7."""
    fixed_snr, fixed_rho = 2.0, 0.7
    fixed_sdx = math.sqrt((1.0 - fixed_rho) / fixed_rho)
    dgp_kwargs = dict(n=cfg["n"], p=cfg["p"], snr=fixed_snr, sd_x=fixed_sdx)

    methods = ["CV_1SE_SVD", "CV_1SE_CCD"]

    print_section_header(f"Part 2/3: lambda2_method Comparison  |  {tag}")
    print(f"  n={cfg['n']}  p={cfg['p']}  s=150  MC={cfg['num_MC']}"
          f"  tFDR={fmt_num(cfg['tFDR'])}  K={cfg['K']}\n"
          f"  fixed operating point: SNR={fmt_num(fixed_snr)}"
          f"  rho={fmt_num(fixed_rho)}\n"
          "  lambda2 methods: CV_1SE_SVD / CV_1SE_CCD"
          "  |  selectors: TENET, TIENET_AUG  (4 runs)\n")

    en_results, ien_results = [], []
    for lm in methods:
        en_results.append(run_gvs_mc_trials(
            "dgp_hastie_snr", dgp_kwargs, cfg["num_MC"], f"TENET+{lm}",
            cfg["tFDR"], cfg["K"], "EN", cfg["corr_max"],
            en_solver="TENET", lambda2_method=lm,
            base_seed=cfg["base_seed"], max_workers=N_WORKERS))
        ien_results.append(run_gvs_mc_trials(
            "dgp_hastie_snr", dgp_kwargs, cfg["num_MC"], f"TIENET_AUG+{lm}",
            cfg["tFDR"], cfg["K"], "IEN", cfg["corr_max"],
            lambda2_method=lm,
            base_seed=cfg["base_seed"], max_workers=N_WORKERS))

    save_mc_1d_method_tables(
        tag, "lambda2_method", "lambda2_method", methods, en_results,
        ien_results, f"SNR={fmt_num(fixed_snr)}, rho={fmt_num(fixed_rho)}",
        cfg, OUT_DIR, f"gvs_{tag}_lambda2_method")


# ==============================================================================
#  Part 3 — hc_linkage comparison
# ==============================================================================

def sim_hastie_part_3(cfg, tag):
    """Compare HAC linkage rules (Single / Complete / Average / WPGMA) for
    TENET and TIENET_AUG at the same fixed operating point as Part 2, with
    lambda2 fixed at CV_1SE_CCD to isolate the linkage effect."""
    fixed_snr, fixed_rho = 2.0, 0.7
    fixed_sdx = math.sqrt((1.0 - fixed_rho) / fixed_rho)
    dgp_kwargs = dict(n=cfg["n"], p=cfg["p"], snr=fixed_snr, sd_x=fixed_sdx)

    linkages = ["Single", "Complete", "Average", "WPGMA"]

    print_section_header(f"Part 3/3: hc_linkage Comparison  |  {tag}")
    print(f"  n={cfg['n']}  p={cfg['p']}  s=150  MC={cfg['num_MC']}"
          f"  tFDR={fmt_num(cfg['tFDR'])}  K={cfg['K']}  lambda2=CV_1SE_CCD\n"
          f"  fixed operating point: SNR={fmt_num(fixed_snr)}"
          f"  rho={fmt_num(fixed_rho)}\n"
          "  linkages: Single / Complete / Average / WPGMA"
          "  |  selectors: TENET, TIENET_AUG  (8 runs)\n")

    en_results, ien_results = [], []
    for lk in linkages:
        en_results.append(run_gvs_mc_trials(
            "dgp_hastie_snr", dgp_kwargs, cfg["num_MC"], f"TENET+{lk}",
            cfg["tFDR"], cfg["K"], "EN", cfg["corr_max"],
            hc_linkage=lk, en_solver="TENET",
            base_seed=cfg["base_seed"], max_workers=N_WORKERS))
        ien_results.append(run_gvs_mc_trials(
            "dgp_hastie_snr", dgp_kwargs, cfg["num_MC"], f"TIENET_AUG+{lk}",
            cfg["tFDR"], cfg["K"], "IEN", cfg["corr_max"],
            hc_linkage=lk,
            base_seed=cfg["base_seed"], max_workers=N_WORKERS))

    save_mc_1d_method_tables(
        tag, "hc_linkage", "hc_linkage", linkages, en_results, ien_results,
        f"SNR={fmt_num(fixed_snr)}, rho={fmt_num(fixed_rho)}",
        cfg, OUT_DIR, f"gvs_{tag}_hc_linkage")


# ==============================================================================
#  main
# ==============================================================================

if __name__ == "__main__":
    print(f"Running with {N_WORKERS} worker processes\n")

    if RUN_PART_1:
        sim_hastie_part_1(CFG, SCENARIO_TAG)
    if RUN_PART_2:
        sim_hastie_part_2(CFG, SCENARIO_TAG)
    if RUN_PART_3:
        sim_hastie_part_3(CFG, SCENARIO_TAG)

    print(f"\n{SCENARIO_TAG} GVS MC simulations complete.")
