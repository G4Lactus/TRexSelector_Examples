# ==============================================================================
# trex_spca_sim_utils.py
# ==============================================================================
#
# Shared simulation utilities for the T-Rex SPCA demo files. Python counterpart
# of cpp/trex_selector_methods/trex_spca/trex_spca_sim_utils.hpp.
#
#   Data generation:
#     generate_sparse_factor_model()  — sparse M-factor model DGP
#                                       (returns dict: X, Z, V, E, actual_snr_db)
#
#   Preprocessing:
#     center_columns()                — column centering (covariance-PCA footing)
#
#   Evaluation:
#     evaluate_spca()                 — PC1 TPR/FDR + cumulative PEV (QR-adjusted)
#
#   MC loop functions (parallel via multiprocessing spawn):
#     run_mc_trials_trex_spca()       — MC loop for T-Rex SPCA
#     run_mc_trials_pca()             — MC loop for the ordinary PCA baseline
#     run_mc_trials_oracle_pca()      — MC loop for oracle thresholded PCA
#
#   Output:
#     save_and_print_spca_mc_results() — table (console + .txt) + tidy .csv
#
# Reproducibility pattern (mirrors cpp):
#   - Each trial's DATA is drawn from a deterministic seed
#     (base_seed_offset + mc*1000, disjoint 1000-wide bands per trial).
#   - The T-Rex selector runs with seed = -1 (hardware entropy per trial), so
#     the computer-generated dummies are independent across trials — required
#     for a valid Monte Carlo FDR estimate.
#
# Support indices are 0-based (Python convention).
# ==============================================================================

import multiprocessing
import os

import numpy as np
import trex_selector_neo.trex_selector_methods as ts
from trex_selector_neo.ml_methods import PCA

N_CORES = min(6, os.cpu_count() or 1)

_SPCA_MODE = {"ActiveSet": ts.SPCAMode.ActiveSet,
              "Thresholded": ts.SPCAMode.Thresholded}
_EN_SOLVER = {"TENET": ts.ENSolverType.TENET,
              "TENET_AUG": ts.ENSolverType.TENET_AUG}


# ==============================================================================
# Data generator (sparse M-factor model)
# ==============================================================================

def generate_sparse_factor_model(n=50, p=100, p1=5, M=3, target_snr=0.0,
                                 overlap_pool_size=30, seed=42):
    """Generate a single dataset from the sparse factor model X = Z V^T + E.

    Z (n x M) has column-wise N(0, sigma_m^2) factors with sigma = {5, 3, 1}
    (falling back to 1 for M > 3). V (p x M) has exactly p1 non-zeros (value
    0.9) per column, drawn without replacement from a shared candidate pool of
    size overlap_pool_size — so different factors' active sets can overlap.
    E is i.i.d. Gaussian noise scaled to match the target SNR in dB.

    Returns dict(X, Z, V, E, actual_snr_db); X is raw (uncentered).
    """
    if not (0 < p1 <= p) or n <= 0:
        raise ValueError("n, p, p1 must be positive with p1 <= p.")
    rng = np.random.default_rng(seed)

    # 1. True factors Z (n x M) with decreasing signal strength
    factor_stds = [5.0, 3.0, 1.0]
    Z = np.column_stack([
        rng.normal(0.0, factor_stds[m] if m < len(factor_stds) else 1.0, n)
        for m in range(M)])

    # 2. Sparse loadings V (p x M): each factor draws p1 indices from the pool
    actual_pool = min(overlap_pool_size, p)
    if p1 > actual_pool:
        raise ValueError("p1 must be <= the overlap pool size.")
    V = np.zeros((p, M))
    for m in range(M):
        V[rng.choice(actual_pool, p1, replace=False), m] = 0.9

    # 3. Signal matrix S = Z V^T
    S = Z @ V.T
    signal_var = S.var(ddof=1)

    # 4. Gaussian noise E scaled to the target SNR
    # SNR_dB = 10 * log10(var_S / var_E)  =>  var_E = var_S / 10^(SNR_dB / 10)
    E = rng.normal(0.0, 1.0, (n, p))
    target_noise_var = signal_var / 10.0 ** (target_snr / 10.0)
    E *= np.sqrt(target_noise_var / E.var(ddof=1))

    X = S + E
    actual_snr_db = 10.0 * np.log10(signal_var / E.var(ddof=1))
    return dict(X=X, Z=Z, V=V, E=E, actual_snr_db=actual_snr_db)


# ==============================================================================
# Preprocessing
# ==============================================================================

def center_columns(X):
    """Center columns (mean subtraction only, NO scaling).

    Covariance-PCA footing shared by every method, matching the legacy R
    reference pipeline and TRexSPCA's internal convention. Column scales carry
    the factor amplitude signal and must not be normalized (z-scoring destroys
    it — see cpp validation_trex_spca_06_handrolled_comparison).
    """
    return X - X.mean(axis=0)


# ==============================================================================
# Single-trial evaluator
# ==============================================================================

def evaluate_spca(X, V_est, Z_est, V_true):
    """Compute per-trial sparse PCA metrics from estimated loadings and scores.

    TPR / FDR are evaluated strictly on the FIRST principal component (PC1):
    PCA's orthogonality constraint mixes the loading supports of PCs 2+ across
    the true factors, making per-component support recovery metrics ambiguous
    beyond the first PC. Cumulative PEV is computed via QR decomposition of
    Z_est across all M components.

    Returns dict(tpr, fdr, pev).
    """
    n = X.shape[0]
    M = V_est.shape[1]

    # 1. TPR / FDR strictly on PC1
    true_supp_pc1 = set(np.flatnonzero(V_true[:, 0] != 0.0))
    est_supp_pc1 = set(np.flatnonzero(np.abs(V_est[:, 0]) > 1e-12))

    tp = len(true_supp_pc1 & est_supp_pc1)
    tpr = 1.0 if not true_supp_pc1 else tp / len(true_supp_pc1)
    fdr = 0.0 if not est_supp_pc1 else (len(est_supp_pc1) - tp) / len(est_supp_pc1)

    # 2. Cumulative PEV across ALL components (QR of Z_est)
    R = np.linalg.qr(Z_est, mode="r")
    total_var = (X ** 2).sum() / (n - 1)
    cum_ev = (np.diag(R)[:M] ** 2).sum() / (n - 1)
    pev = cum_ev / total_var if total_var > 0.0 else 0.0

    return dict(tpr=tpr, fdr=fdr, pev=pev)


# ==============================================================================
# Per-trial workers (top-level, importable by multiprocessing spawn workers)
# ==============================================================================

def _trial_trex_spca(task):
    seed, dgp_kwargs, tFDR, mode, lambda_2, en_solver = task
    dat = generate_sparse_factor_model(seed=seed, **dgp_kwargs)
    X = center_columns(dat["X"])
    M = dat["V"].shape[1]

    ctrl = ts.TRexSPCAControlParameter()
    ctrl.mode = _SPCA_MODE[mode]
    ctrl.en_solver = _EN_SOLVER[en_solver]
    ctrl.gvs_ctrl.lambda2_method = ts.LambdaSelectionMethod.CV_1SE_CCD
    ctrl.gvs_ctrl.lambda_2 = lambda_2  # < 0: auto via CV; 0: no ridge; > 0: fixed

    # Selector seed -1: hardware entropy per trial (dummies independent across
    # trials). The selector maps X zero-copy — pass a Fortran-ordered copy.
    sel = ts.TRexSPCASelector(np.asfortranarray(X), M, tFDR, ctrl, -1, False)
    sel.select()
    res = sel.getResult()

    m = evaluate_spca(X, np.asarray(res.V), np.asarray(res.Z), dat["V"])
    return m["tpr"], m["fdr"], m["pev"]


def _trial_pca(task):
    seed, dgp_kwargs = task
    dat = generate_sparse_factor_model(seed=seed, **dgp_kwargs)
    X = center_columns(dat["X"])
    M = dat["V"].shape[1]

    # X is already centered -> center=False, normalize=False (covariance PCA).
    pca = PCA(center=False, normalize=False)
    res = pca.fit(np.asfortranarray(X), M)

    m = evaluate_spca(X, np.asarray(res.V), np.asarray(res.Z), dat["V"])
    return m["tpr"], m["fdr"], m["pev"]


def _trial_oracle_pca(task):
    seed, dgp_kwargs, p1 = task
    dat = generate_sparse_factor_model(seed=seed, **dgp_kwargs)
    X = center_columns(dat["X"])
    p = X.shape[1]
    M = dat["V"].shape[1]

    pca = PCA(center=False, normalize=False)
    ord_res = pca.fit(np.asfortranarray(X), M)
    V_ord = np.asarray(ord_res.V)

    # Keep top p1 entries per loading vector, then normalize
    V_oracle = np.zeros((p, M))
    for comp in range(M):
        v_m = V_ord[:, comp]
        thresh = np.sort(np.abs(v_m))[p - p1]
        keep = np.abs(v_m) >= thresh
        V_oracle[keep, comp] = v_m[keep]
        V_oracle[:, comp] /= np.linalg.norm(V_oracle[:, comp])
    Z_oracle = X @ V_oracle

    m = evaluate_spca(X, V_oracle, Z_oracle, dat["V"])
    return m["tpr"], m["fdr"], m["pev"]


# ==============================================================================
# MC aggregation helper
# ==============================================================================

def _run_mc(worker, tasks, label, n_cores):
    print(f"  {label} — Running {len(tasks)} MC trials ...", flush=True)
    if n_cores > 1:
        with multiprocessing.get_context("spawn").Pool(n_cores) as pool:
            res = pool.map(worker, tasks)
    else:
        res = [worker(t) for t in tasks]

    tpr = np.array([r[0] for r in res])
    fdr = np.array([r[1] for r in res])
    pev = np.array([r[2] for r in res])
    out = dict(avg_tpr=tpr.mean(), avg_fdr=fdr.mean(), avg_pev=pev.mean(),
               sd_tpr=tpr.std(ddof=1) if len(tpr) > 1 else 0.0,
               sd_fdr=fdr.std(ddof=1) if len(fdr) > 1 else 0.0)
    print(f"  {label} — done. TPR={out['avg_tpr']:.3f}  FDR={out['avg_fdr']:.3f}\n",
          flush=True)
    return out


# ==============================================================================
# Parallel MC inner loops
# ==============================================================================

def run_mc_trials_trex_spca(num_MC, progress_label, dgp_kwargs, tFDR, mode,
                            base_seed_offset, lambda_2=-1.0, en_solver="TENET",
                            n_cores=N_CORES):
    """Run num_MC parallel T-Rex SPCA trials.

    Per trial: the DGP draws data from seed base_seed_offset + mc*1000; the
    pipeline centers X once (covariance PCA); this same X is shared by the
    selector and by evaluate_spca(). mode is "ActiveSet" or "Thresholded";
    en_solver is "TENET" or "TENET_AUG". Returns dict(avg_tpr, avg_fdr,
    avg_pev, sd_tpr, sd_fdr).
    """
    tasks = [(base_seed_offset + mc * 1000, dgp_kwargs, tFDR, mode,
              lambda_2, en_solver) for mc in range(num_MC)]
    return _run_mc(_trial_trex_spca, tasks, progress_label, n_cores)


def run_mc_trials_pca(num_MC, progress_label, dgp_kwargs, base_seed_offset,
                      n_cores=N_CORES):
    """Run num_MC parallel ordinary PCA trials.

    Ordinary PCA retains all p loadings per component (no sparsity), so
    FDR ~ (p - p1) / p and TPR = 1 regardless of SNR; the metric of interest
    is the cumulative PEV.
    """
    tasks = [(base_seed_offset + mc * 1000, dgp_kwargs) for mc in range(num_MC)]
    return _run_mc(_trial_pca, tasks, progress_label, n_cores)


def run_mc_trials_oracle_pca(num_MC, progress_label, dgp_kwargs, p1,
                             base_seed_offset, n_cores=N_CORES):
    """Run num_MC parallel oracle-thresholded PCA trials.

    Oracle PCA keeps exactly the top p1 entries (by absolute value) of each
    ordinary-PCA loading vector — it knows the true support cardinality. This
    bounds what any data-driven thresholding method can achieve.
    """
    tasks = [(base_seed_offset + mc * 1000, dgp_kwargs, p1)
             for mc in range(num_MC)]
    return _run_mc(_trial_oracle_pca, tasks, progress_label, n_cores)


# ==============================================================================
# Output
# ==============================================================================

def save_and_print_spca_mc_results(out_dir, file_stem, num_MC, snr_values,
                                   method_names, fdr_map, tpr_map, pev_map):
    """Print MC-averaged results as an aligned table (console + .txt) and
    write a tidy long-format CSV (method, metric, snr_db, value).

    fdr_map / tpr_map / pev_map are dicts: method name -> list of values over
    the SNR grid.
    """
    os.makedirs(out_dir, exist_ok=True)

    lines = []
    lines.append("")
    lines.append("=" * 70)
    lines.append(f"=== T-Rex SPCA Results (averaged over {num_MC} MC trials) ===")
    lines.append("=" * 70)
    lines.append("")

    hdr = f"{'Method':<22}{'Metric':<8}{'SNR(dB)':>8}"
    hdr += "".join(f"{snr:>10.1f}" for snr in snr_values)
    lines.append(hdr)
    lines.append("-" * (38 + 10 * len(snr_values)))

    for name in method_names:
        for metric, mmap in (("FDR", fdr_map), ("TPR", tpr_map), ("PEV", pev_map)):
            first = metric == "FDR"
            row = f"{name if first else '':<22}{metric:<8}{'':>8}"
            row += "".join(f"{v:>10.4f}" for v in mmap[name])
            lines.append(row)
        lines.append("")

    text = "\n".join(lines) + "\n"
    print(text, end="", flush=True)
    txt_path = os.path.join(out_dir, file_stem + ".txt")
    with open(txt_path, "w") as f:
        f.write(text)

    csv_path = os.path.join(out_dir, file_stem + ".csv")
    with open(csv_path, "w") as f:
        f.write("method,metric,snr_db,value\n")
        for name in method_names:
            for i, snr in enumerate(snr_values):
                f.write(f"{name},FDR,{snr:.6f},{fdr_map[name][i]:.6f}\n")
                f.write(f"{name},TPR,{snr:.6f},{tpr_map[name][i]:.6f}\n")
                f.write(f"{name},PEV,{snr:.6f},{pev_map[name][i]:.6f}\n")
    print(f"[Info] CSV results saved to:             {csv_path}")
    print(f"[Info] Results successfully saved to: {txt_path}\n")
