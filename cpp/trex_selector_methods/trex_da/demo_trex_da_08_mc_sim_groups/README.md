# Demo 08: DA-TRex on Multi-Level Prior-Groups Data

## Purpose

Test DA-TRex on the structurally richest dependency model in this folder: a **three-level hierarchical latent-factor structure** with a non-exchangeable Toeplitz leaf layer, using explicit prior group labels rather than an automatically-discovered dendrogram cut.

---

## Data Generation Parameters (`dgp_groups_toeplitz_leaf`)

$$
X_{i,j} = \sum_{x=1}^{L-1} \sqrt{\rho_x - \rho_{x+1}}\, f_{g_x(j)} + \sqrt{\rho_0 - \rho_1}\, u_{g_0(j),j} + \sqrt{1-\rho_0}\,\varepsilon_j,
$$

where coarser levels $x\ge1$ use standard latent factors, and the finest level ($x=0$) uses a **Toeplitz-correlated leaf block** ($\mathrm{Cov}(u_r,u_s)=\phi^{|r-s|}$) instead of a single shared factor.

- **Three nested grouping levels**: Level 1 = groups of $10$ ($p/10$ groups), Level 2 = groups of $50$ ($p/50$ groups), Level 3 = groups of $250$ ($p/250$ groups).
- Cumulative correlation levels (fine → coarse): $\rho=\{0.55, 0.25, 0.10\}$; Toeplitz leaf decay $\phi_{\text{leaf}}=0.50$.
- $n=300$, $p=1000$, $s=10$, `Random` support, `max_gap=20` (unused directly by this DGP but present in the shared config).
- Prior groups are passed directly to the selector via `da_ctrl.prior_groups` (the 3-level `groups` structure) and `da_ctrl.rho_grid_labels` (the $\rho$ levels) — this demo supplies the hierarchy explicitly rather than having DA-TRex discover it via HAC clustering.

---

## Control Parameters

```
K = 20
tFDR = 0.2
MC = 200
base_seed = 2026
SNR grid = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0}
```

Solvers: TLARS, TAFS, TOMP.

---

## Output Files

Written to `simulation_results/`:

- `da_trex_mc_da_groups_toeplitz_leaf.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_da/demo_trex_da_08_mc_sim_groups/demo_trex_da_08_mc_sim_groups
```

---

## Real Results (excerpt: SNR sweep, TLARS and TAFS)

| SNR | TLARS FDR | TLARS TPR | TAFS FDR | TAFS TPR |
|---|---|---|---|---|
| 0.10 | 0.0742 | 0.0070 | 0.0829 | 0.0090 |
| 0.50 | 0.3540 | 0.3685 | 0.2106 | 0.3340 |
| 1.00 | 0.3521 | 0.6855 | 0.1138 | 0.7305 |
| 2.00 | 0.3488 | 0.7990 | 0.0769 | 0.9465 |
| 5.00 | 0.3314 | 0.8000 | 0.1048 | 0.9400 |

---

## Interpretation

- TLARS's FDR climbs well **above** the $\mathrm{tFDR}=0.2$ target from $\mathrm{SNR}=0.5$ onward and plateaus around $0.33$–$0.35$ — the highest sustained FDR overshoot of any real-data demo in this folder — suggesting that this three-level nested dependency structure is genuinely harder for DA-TRex to control than the single-level structures in Demos 01–07, at least with the TLARS solver.
- **TAFS notably outperforms TLARS here**: at $\mathrm{SNR}=2.0$, TAFS achieves FDR $0.077$ (well controlled) with TPR $0.9465$ (higher than TLARS's $0.799$) — a rare case in this project where solver choice materially changes both FDR control and power, rather than the usual "solvers agree, only speed differs" pattern seen in the classical T-Rex demos. This is worth investigating further if you extend this demo.
- Because prior groups are supplied directly (not HAC-discovered), any residual FDR inflation here reflects genuine difficulty in the underlying multi-level dependency-aware correction itself, not a clustering-quality artifact — this differentiates Demo 08 from the `BT`-based Demos 04–07, where clustering quality is an additional confound.
- This demo currently only exercises the SNR sweep; extending it with $\rho$-level or group-size sweeps (analogous to Demos 04–07's $Q$/$M$ sweeps) would help clarify which aspect of the nested structure drives the elevated FDR.

---

**Last updated**: 2026-07-08
