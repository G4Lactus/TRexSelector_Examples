# Demo 07: DA-TRex (BT) on Heavy-Tailed Block + White-Noise Data

## Purpose

Combine Demo 05's white-noise dilution with Demo 06's heavy-tailed robustness study: heavy-tailed AR(1)-Toeplitz blocks embedded in a larger design padded with i.i.d. Student-$t$ white-noise columns. Corresponds to R reference `demo_trex_da_07_bt_ht_block_white_sweeps.R`.

---

## Data Generation Parameters (`dgp_ht_block_white`)

- $M$ blocks of size $Q$ with heavy-tailed ($t(\nu)$) AR(1)-Toeplitz within-block covariance (always heavy-tailed for $X$, unlike Demo 06 where Gaussian $X$ is also an option), plus $p_{\text{white}}$ appended i.i.d. $t(\nu)$ white-noise columns; **total dimension fixed at $p_{\text{total}}=500$**.
- Response noise is Gaussian or $t(\nu)$-distributed depending on the scenario (same `s1_Gauss`/`s2_Heavy` labeling as Demo 06).
- Active variables ($s=M$, one per AR block) lie in the heavy-tailed AR part only.
- Base parameters: $n=150$, $p_{\text{total}}=500$, $M=5$, $Q=5$ ($p_{ar}=25$, $p_{\text{white}}=475$, $s=5$), $\rho=0.8$, amplitude $1.0$, $\nu=3.0$, $\mathrm{tFDR}=0.2$, $K=20$, $\mathrm{MC}=200$, seed $2026$.

---

## Five Sections (each looped over noise scenario × linkage)

1. **SNR sweep** — $\mathrm{SNR}\in\{0.1,0.2,0.5,0.6,1.0,2.0,5.0\}$.
2. **$\rho$ sweep** — $\rho\in\{0.0,0.1,\dots,0.9\}$.
3. **$Q$ sweep** — $Q\in\{5,\dots,50\}$.
4. **$M$ sweep** — $M\in\{1,\dots,30\}$.
5. **$\mathrm{tFDR}$ sweep** — $\mathrm{tFDR}\in\{0.05,\dots,0.50\}$.

(No dedicated linkage-only section here, unlike Demo 06 — linkage is still swept as an outer loop within each of the five sections above.)

---

## Output Files

Written to `simulation_results/` (48 files): naming pattern `da_trex_mc_da_ht_white_{snr|rho|Q|M|tFDR}_{s1_Gauss|s2_Heavy}_{Single|Complete|Average}.{txt,csv}`.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_da_07_mc_sim_bt_ht_block_white
```

---

## Real Results (excerpt: SNR sweep, Heavy noise, Single linkage, TLARS)

| SNR | FDR | TPR |
|---|---|---|
| 0.10 | 0.1125 | 0.0040 |
| 0.50 | 0.0957 | 0.2980 |
| 1.00 | 0.0680 | 0.5950 |
| 2.00 | 0.0449 | 0.7720 |
| 5.00 | 0.0367 | 0.8040 |

---

## Interpretation

- Unlike Demo 06 (heavy-tailed blocks *without* white-noise dilution), FDR here **decreases** as SNR increases and ends up quite low ($\approx0.037$–$0.045$) at high SNR — closer to the white-noise-diluted pattern of Demo 05 than to Demo 06's roughly flat, higher FDR — suggesting the large ($475$-column) white-noise padding helps FDR control even under heavy tails, similar to how it helped in the Gaussian case (Demo 05 vs. Demo 04).
- TPR reaches $\approx0.80$ by $\mathrm{SNR}=5.0$, notably higher than Demo 06's $\approx0.63$ plateau at the same nominal SNR range — again consistent with the white-noise dilution making the true AR(1) blocks comparatively easier to isolate than in Demo 06's smaller, all-block design.
- Read this demo as the "combined effects" check: compare against Demo 04 (Gaussian, no white noise), Demo 05 (Gaussian, with white noise), and Demo 06 (heavy-tailed, no white noise) to decompose how much each factor — heavy tails vs. white-noise dilution — independently affects FDR/TPR, and whether their effects on this data compound or partially offset.
- As in Demo 06, the $\mathrm{tFDR}$ sweep (Section 5) is useful for checking whether realized FDR tracks the target level consistently across this different (white-noise-diluted) design.

---

**Last updated**: 2026-07-04
