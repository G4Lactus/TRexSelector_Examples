# ==============================================================================
# demo_trex_spca_01_mc_sim.py
# ==============================================================================
#
# Monte Carlo simulation comparing T-Rex SPCA against PCA baselines.
# Python counterpart of cpp demo_trex_spca_01_mc_sim.cpp.
#
# DGP: sparse M-factor model — X = Z * V^T + E.
#      n=50, p=100, M=3 latent factors, p1=5 active loadings per factor,
#      overlap_pool_size=30 (shared candidate-index pool across factors).
#
# Section 1: SNR sweep over the snr_values grid set in main()
#            (currently the single point -10 dB, num_MC=200).
#   Methods compared:
#     1. OrdPCA             — ordinary PCA (baseline; no sparsity constraint)
#     2. OraclePCA          — thresholded PCA with known support cardinality p1
#     3. TRexSPCA-EN-Act    — T-Rex+GVS(EN/TENET, Gram-based) + ActiveSet
#     4. TRexSPCA-ENaug-Act — T-Rex+GVS(EN/TENETAug, augmented) + ActiveSet
#     5. TRexSPCA-EN-Thr    — T-Rex+GVS(EN/TENET, Gram-based) + Thresholded
#     6. TRexSPCA-ENaug-Thr — T-Rex+GVS(EN/TENETAug, augmented) + Thresholded
#
#   All methods see the same center-only X (covariance-PCA footing, as in the
#   legacy R reference); the per-PC selector uses the default L2 scaling.
#
# Metrics (num_MC-averaged; TPR/FDR on PC1 support only — see evaluate_spca):
#   - FDR  — false discovery rate of PC1's estimated loading support
#   - TPR  — true positive rate  of PC1's estimated loading support
#   - PEV  — cumulative percentage of explained variance (all M components)
#
# Support indices are 0-based (Python convention).
#
# Usage: python demo_trex_spca_01_mc_sim.py [num_worker_cores]
# ==============================================================================

import sys
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from trex_spca_sim_utils import (N_CORES, run_mc_trials_pca,
                                 run_mc_trials_oracle_pca,
                                 run_mc_trials_trex_spca,
                                 save_and_print_spca_mc_results)

N_WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else N_CORES


# ==============================================================================
# Monte Carlo Simulation — SNR sweep
# ==============================================================================

def demo_trex_spca_mc_snr_sweep(cfg, snr_values, stem):
    """Run the SNR-sweep MC simulation and save results.

    cfg holds the simulation parameters (n, p, p1, M, overlap_pool_size,
    tFDR, lambda_2, num_MC, base_seed); snr_values is the SNR grid in dB;
    stem is the output file base name (no folder, no extension).
    """
    print("================================================================")
    print("  T-Rex SPCA MC Simulation — SNR Sweep")
    print("================================================================")
    print(f"  n={cfg['n']}  p={cfg['p']}  M={cfg['M']}  p1={cfg['p1']}"
          f"  tFDR={cfg['tFDR']}  num_MC={cfg['num_MC']}"
          f"  ({N_WORKERS} workers)\n")

    method_names = ["OrdPCA", "OraclePCA",
                    "TRexSPCA-EN-Act", "TRexSPCA-ENaug-Act",
                    "TRexSPCA-EN-Thr", "TRexSPCA-ENaug-Thr"]

    fdr_map = {name: [] for name in method_names}
    tpr_map = {name: [] for name in method_names}
    pev_map = {name: [] for name in method_names}

    def store(name, r):
        fdr_map[name].append(r["avg_fdr"])
        tpr_map[name].append(r["avg_tpr"])
        pev_map[name].append(r["avg_pev"])

    dgp_base = dict(n=cfg["n"], p=cfg["p"], p1=cfg["p1"], M=cfg["M"],
                    overlap_pool_size=cfg["overlap_pool_size"])

    for snr in snr_values:
        offset = cfg["base_seed"] + int(round(snr))
        snr_lbl = f"SNR={snr:.1f} dB"
        dgp_kwargs = dict(dgp_base, target_snr=snr)

        print("--------------------------------------------------")
        print(f"{snr_lbl}\n")

        # 1. Ordinary PCA
        store("OrdPCA", run_mc_trials_pca(
            cfg["num_MC"], f"{snr_lbl} [OrdPCA]", dgp_kwargs, offset,
            N_WORKERS))

        # 2. Oracle Thresholded PCA
        store("OraclePCA", run_mc_trials_oracle_pca(
            cfg["num_MC"], f"{snr_lbl} [OraclePCA]", dgp_kwargs, cfg["p1"],
            offset, N_WORKERS))

        # 3. T-Rex SPCA — EN (Gram-based TENET) + ActiveSet
        store("TRexSPCA-EN-Act", run_mc_trials_trex_spca(
            cfg["num_MC"], f"{snr_lbl} [TRexSPCA-EN-Act]", dgp_kwargs,
            cfg["tFDR"], "ActiveSet", offset, cfg["lambda_2"], "TENET",
            N_WORKERS))

        # 4. T-Rex SPCA — ENaug (augmented TENET) + ActiveSet
        store("TRexSPCA-ENaug-Act", run_mc_trials_trex_spca(
            cfg["num_MC"], f"{snr_lbl} [TRexSPCA-ENaug-Act]", dgp_kwargs,
            cfg["tFDR"], "ActiveSet", offset, cfg["lambda_2"], "TENET_AUG",
            N_WORKERS))

        # 5. T-Rex SPCA — EN (Gram-based TENET) + Thresholded
        store("TRexSPCA-EN-Thr", run_mc_trials_trex_spca(
            cfg["num_MC"], f"{snr_lbl} [TRexSPCA-EN-Thr]", dgp_kwargs,
            cfg["tFDR"], "Thresholded", offset, cfg["lambda_2"], "TENET",
            N_WORKERS))

        # 6. T-Rex SPCA — ENaug (augmented TENET) + Thresholded
        store("TRexSPCA-ENaug-Thr", run_mc_trials_trex_spca(
            cfg["num_MC"], f"{snr_lbl} [TRexSPCA-ENaug-Thr]", dgp_kwargs,
            cfg["tFDR"], "Thresholded", offset, cfg["lambda_2"], "TENET_AUG",
            N_WORKERS))

    save_and_print_spca_mc_results(
        out_dir=os.path.join(_THIS_DIR, "simulation_results", "data"),
        file_stem=stem, num_MC=cfg["num_MC"],
        snr_values=snr_values, method_names=method_names,
        fdr_map=fdr_map, tpr_map=tpr_map, pev_map=pev_map)


# ==============================================================================
# Main
# ==============================================================================

def main():
    # =========================================================
    # Section 1: SNR sweep
    # =========================================================
    if True:
        cfg = dict(
            n=50,
            p=100,
            p1=5,
            M=3,
            overlap_pool_size=30,
            tFDR=0.10,
            num_MC=200,
            base_seed=42,
            lambda_2=-1.0,  # -1 triggers k-fold CV for lambda_2 selection
        )

        snr_values = [-10.0]  # TEMP fast-validation

        demo_trex_spca_mc_snr_sweep(cfg, snr_values, "demo_trex_spca_01_mc_sim")


if __name__ == "__main__":
    main()
