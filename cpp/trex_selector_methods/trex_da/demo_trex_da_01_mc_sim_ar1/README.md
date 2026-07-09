# Demo 01: DA-TRex on AR(1) Toeplitz Data

## Purpose

The foundational DA-TRex scenario: recover a sparse support under AR(1) Toeplitz column correlation, comparing **DA-TRex (`AR1`)** directly against the **classical (no-DA) T-Rex selector** on the same data, across an SNR sweep, a $\rho$ sweep, and a 2D gap-vs-$\rho$ study of when the correction breaks down.

---

## Data Generation Parameters (`dgp_ar1`)

$$
X_{i,j} = \rho\, X_{i,j-1} + \sqrt{1-\rho^2}\, \eta_{i,j}, \qquad \Sigma_{jk} = \rho^{|j-k|}.
$$

- $n=300$, $p=1000$, $s=10$ active variables, coefficient amplitude $3.0$.
- $\mathrm{tFDR}=0.2$.

---

## Three Parts

### Part 1 — SNR sweep
$\rho=0.7$ fixed; $\mathrm{SNR}\in\{0.1,0.2,0.5,1.0,2.0,5.0\}$; two support placements: `CappedSpread(max_gap=20)` and `Random` (redrawn per trial).

### Part 2 — $\rho$ sweep
$\mathrm{SNR}=2.0$ fixed; $\rho\in\{0.0,0.1,\dots,0.9\}$; same two support placements.

### Part 3 — 2D gap × $\rho$ sweep
Explores the **$\kappa$-boundary**: the DA+AR1 correction window half-width $\kappa = \lceil \log(0.02)/\log(\rho) \rceil$. When the `CappedSpread` gap between active variables falls below $\kappa$, active variables fall inside each other's correction windows and TPR is expected to collapse.
- `gap_grid = {100, 50, 20, 15, 10, 5, 1}`, `rho_grid = {0.0, 0.1, ..., 0.9}`, $\mathrm{SNR}=2.0$ fixed.

All parts run $K=20$, $\mathrm{MC}=200$ per grid point, solvers TLARS/TAFS/TOMP, plus a **base T-Rex (no DA)** comparison row.

---

## Output Files

Written to `simulation_results/`:

- `da_trex_mc_da_ar1_snr_capped.txt` / `.csv`
- `da_trex_mc_da_ar1_snr_random.txt` / `.csv`
- `da_trex_mc_da_ar1_rho_capped.txt` / `.csv`
- `da_trex_mc_da_ar1_rho_random.txt` / `.csv`
- `da_trex_mc_da_ar1_gap_rho.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1
```

---

## Real Results (excerpt: SNR sweep, CappedSpread support, TLARS)

| SNR | DA-TRex FDR | DA-TRex TPR | Base T-Rex FDR | Base T-Rex TPR |
|---|---|---|---|---|
| 0.10 | 0.0175 | 0.0040 | 0.0567 | 0.0090 |
| 0.50 | 0.1163 | 0.2555 | 0.2804 | 0.4270 |
| 2.00 | 0.0593 | 0.8070 | 0.2899 | 0.9965 |
| 5.00 | 0.0748 | 0.8490 | 0.2856 | 1.0000 |

---

## Interpretation

- **DA correction substantially reduces FDR** relative to base T-Rex at every SNR level shown — at $\mathrm{SNR}=2.0$, DA-TRex's FDR ($0.059$) is roughly a fifth of base T-Rex's ($0.29$) — but this comes at a real TPR cost at higher SNR ($0.807$ vs. $0.9965$), since the DA correction actively excludes some correlated neighbors of true positives that base T-Rex would happily (over-)select.
- At very low SNR ($0.1$), both methods have low TPR since there is little signal to detect either way, but DA-TRex still keeps FDR much lower.
- The 2D gap×$\rho$ study is the key diagnostic for understanding *why* DA-TRex sometimes loses TPR: once active variables are spaced closer together than the $\kappa$-boundary implied by $\rho$, they start suppressing each other under the DA correction — this is an expected, structural property of the method, not a bug.
- Compare `CappedSpread` vs. `Random` support results to see whether deterministic evenly-spaced supports behave differently from randomly placed ones at the same nominal spacing.

---

**Last updated**: 2026-07-08
