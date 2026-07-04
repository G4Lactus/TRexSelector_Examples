# Demo 01: T-Rex+GVS on Hastie Equicorrelated Blocks

## Purpose

Evaluate grouped-signal recovery under strong equicorrelation, using the classical three-group design of Zou & Hastie (2005). This is the canonical starting point for the T-Rex+GVS demos: all variables within an active group are themselves active, so the question is purely about **group-level recovery**, not within-group sparse support identification.

---

## Data Generation Parameters (`make_hastie_dgp`)

$$
X_{ij} = Z_{i,\,g(j)} + \sigma_x\, \xi_{ij}, \qquad \xi_{ij} \sim \mathcal{N}(0,1).
$$

- Three latent factors $Z_{\cdot,k} \sim \mathcal{N}(0,1)$ for $k \in \{0,1,2\}$.
- **Group 0** (columns 0–49), **Group 1** (columns 50–99), **Group 2** (columns 100–149) — all active, $\beta_j = 3$.
- **Background** (columns 150–499): i.i.d. $\mathcal{N}(0,1)$, each a singleton group.
- Within-group correlation $\rho = 1/(1+\sigma_x^2)$; the default $\sigma_x = \sqrt{0.01}$ gives $\rho \approx 0.99$.
- $n = 200$, $p = 500$, $s = 150$.
- SNR calibrated via $\sigma_y = \sqrt{\mathrm{Var}(X\beta)/\mathrm{SNR}}$.

---

## Control Parameters

```
K = 20                      # Dummy experiments per T-loop iteration
tFDR = 0.1                  # Target FDR
corr_max = 0.98             # HAC auto-clustering correlation threshold
hc_linkage = Single          # Single-linkage HAC
lambda2_method = CV_1SE_CCD # Elastic-net penalty selection
MC = 200                    # Monte Carlo repetitions per grid point
```

---

## Methods Compared

- **EN** (`GVSType::EN`, `ENSolverType::TENET`)
- **EN+AUG** (`GVSType::EN`, `ENSolverType::TENET_AUG`)
- **IEN** (`GVSType::IEN`)

---

## Three Parts

1. **SNR sweep** at fixed $\sigma_x = \sqrt{0.01}$ ($\rho \approx 0.99$): $\mathrm{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$.
2. **$\rho$ sweep** at fixed $\mathrm{SNR} = 2.0$: $\rho \in \{0.10, 0.20, \dots, 0.99\}$ (derived via $\sigma_x = \sqrt{(1-\rho)/\rho}$).
3. **2-D SNR $\times$ $\rho$ grid**: 5 SNR levels $\times$ 6 $\rho$ levels.

---

## Output Files

Written to `simulation_results/`:

- `gvs_Hastie_snr.txt` / `.csv` — SNR-sweep mean FDP/TPP table (Part 1).
- `gvs_Hastie_rho.txt` / `.csv` — $\rho$-sweep mean FDP/TPP table (Part 2).
- `gvs_Hastie_2d.txt` / `.csv` — 2-D SNR $\times$ $\rho$ grid (Part 3).
- `gvs_Hastie_hc_linkage.txt` / `.csv` — supplementary comparison across HAC linkage methods.
- `gvs_Hastie_lambda2_method.txt` / `.csv` — supplementary comparison across $\lambda_2$ selection methods.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_gvs_01_mc_sim_hastie_en_blocks
```

---

## Interpretation

- FDR is expected to stay near the $\mathrm{tFDR}=0.1$ target across the SNR grid for **EN** and **EN+AUG**, with TPR rising from weak recovery at low SNR toward full recovery ($\mathrm{TPR}\approx 1$) once $\mathrm{SNR} \gtrsim 1$.
- **IEN** is expected to behave much more conservatively on this design (very low FDR, low TPR that grows only slowly with SNR) compared to the EN-based methods — this conservative pattern is not specific to this scenario and recurs across the other block demos in this folder as well.
- The $\rho$ sweep is expected to show stable performance across most within-group correlation levels, since all three methods use HAC-based grouping to aggregate evidence across strongly correlated columns.
- The supplementary linkage and $\lambda_2$-method comparisons are useful for checking that the qualitative conclusions above are not sensitive to those implementation choices.

---

**Last updated**: 2026-07-04
