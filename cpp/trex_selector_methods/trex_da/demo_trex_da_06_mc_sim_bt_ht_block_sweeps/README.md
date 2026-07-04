# Demo 06: DA-TRex (BT) on Heavy-Tailed Block Data

## Purpose

Stress-test `BT`-based DA-TRex under **heavy-tailed (Student-t)** predictors and/or noise, using a block-diagonal Toeplitz design. Studies both statistical robustness (Gaussian vs. heavy-tailed scenarios) and sensitivity across five sweep dimensions plus a dedicated linkage-method sweep. Corresponds to R reference `demo_trex_da_06_bt_heavy_tailed_block_sweeps.R`.

---

## Data Generation Parameters (`dgp_block_toeplitz_hvt`)

- $G=M$ blocks of size $g=Q$, each with within-block Toeplitz covariance $\Sigma_{ij}=\rho^{|i-j|}$ (Cholesky-generated); between blocks, zero correlation.
- Rows optionally follow a multivariate Student-$t(\nu)$ distribution (heavy-tailed $X$); response noise is either Gaussian or $t(\nu)$-distributed depending on the scenario.
- Support is one active variable per block (`OnePerBlock`, internal — no separate support/amplitude arguments).
- Base parameters: $n=150$, $M=5$, $Q=5$ ($p=25$, $s=5$), $\rho=0.8$, $\nu=3.0$, $\mathrm{tFDR}=0.2$, $K=20$, $\mathrm{MC}=200$, seed $2026$.

---

## Six Sections (each looped over noise scenario × linkage)

Outer loops: noise scenario $\in$ {**s1_Gauss** (Gaussian noise), **s2_Heavy** ($t(\nu{=}3)$ noise)} $\times$ linkage $\in$ {Single, Complete, Average}.

1. **SNR sweep** — $\mathrm{SNR}\in\{0.1,0.2,0.5,0.6,1.0,2.0,5.0\}$.
2. **$\rho$ sweep** — $\rho\in\{0.0,0.1,\dots,0.9\}$.
3. **$Q$ sweep** — $Q\in\{5,\dots,50\}$.
4. **$M$ sweep** — $M\in\{1,\dots,30\}$.
5. **$\mathrm{tFDR}$ sweep** — $\mathrm{tFDR}\in\{0.05,\dots,0.50\}$.
6. **Linkage sweep** — dedicated comparison of Single/Complete/Average at fixed base parameters.

---

## Output Files

Written to `simulation_results/` (60 files — the largest output set of any demo in this project): naming pattern `da_trex_mc_da_ht_{snr|rho|Q|M|tFDR|linkage}_{s1_Gauss|s2_Heavy}_{Single|Complete|Average}.{txt,csv}` (the dedicated linkage-sweep section omits the trailing linkage suffix since linkage itself is the swept variable: `da_trex_mc_da_ht_linkage_{s1_Gauss|s2_Heavy}.{txt,csv}`).

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_da_06_mc_sim_bt_ht_block_sweeps
```

---

## Real Results (excerpt: SNR sweep, Heavy noise, Single linkage, TLARS)

| SNR | FDR | TPR |
|---|---|---|
| 0.10 | 0.0729 | 0.0240 |
| 0.50 | 0.1917 | 0.4240 |
| 1.00 | 0.1961 | 0.6130 |
| 2.00 | 0.1840 | 0.6250 |
| 5.00 | 0.2052 | 0.6260 |

---

## Interpretation

- FDR hovers **around or somewhat above** the $\mathrm{tFDR}=0.2$ target across most of the SNR range under heavy-tailed noise — noticeably less tightly controlled than the Gaussian block demos (04–05), reflecting the added estimation difficulty from occasional large outliers in $t(3)$-distributed data.
- TPR plateaus around $0.61$–$0.63$ from $\mathrm{SNR}=1.0$ upward — lower than Demo 04's Gaussian block-AR(1) plateau ($\approx0.76$–$0.79$) at comparable design size, consistent with heavy tails making the true signal harder to separate from noise even at high nominal SNR.
- Compare the **s1_Gauss** vs. **s2_Heavy** scenario files directly (same base parameters, only the noise/predictor distribution differs) to isolate exactly how much heavy tails cost in FDR control and power.
- The dedicated **linkage sweep** (Section 6) is worth checking on its own: with heavy-tailed data, correlation estimates used by HAC clustering are noisier, so linkage-method choice may matter more here than in the cleaner Gaussian block demos.
- The $\mathrm{tFDR}$ sweep (Section 5) is unique to this demo (and Demo 07) — use it to check whether the realized FDR tracks the target linearly across a range of target levels, or whether the heavy-tailed setting causes systematic over- or under-shooting at particular target values.

---

**Last updated**: 2026-07-04
