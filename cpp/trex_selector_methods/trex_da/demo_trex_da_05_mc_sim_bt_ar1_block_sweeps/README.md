# Demo 05: DA-TRex (BT) on Block AR(1) + White-Noise Data

## Purpose

Extend Demo 04's block-diagonal AR(1) design by appending a large block of **i.i.d. white-noise columns**, testing whether `BT`-based DA-TRex can still isolate the true AR(1) blocks when they are a small fraction of a much larger, otherwise-uncorrelated design. Corresponds to R reference `demo_trex_da_06_bt_ar1_plus_white_block_sweeps.R` (numbered "06" in the R suite but "05" in this C++ folder — a naming lineage quirk, not a bug).

---

## Data Generation Parameters (`dgp_ar1_block_white`)

- $M$ AR(1) blocks of size $Q$ (as in Demo 04, $p_{ar} = M \times Q$), plus $p_{\text{white}}$ appended i.i.d. $\mathcal{N}(0,1)$ columns; **total dimension fixed at $p_{\text{total}} = p_{ar} + p_{\text{white}} = 1000$**.
- Active variables ($s=M$, one per AR block) lie only in the AR part.
- Base parameters: $n=300$, $p_{\text{total}}=1000$, $M=5$, $Q=5$ ($p_{ar}=25$, $p_{\text{white}}=975$, $s=5$), $\rho=0.7$, amplitude $1.0$, $\mathrm{tFDR}=0.1$, $K=20$, $\mathrm{MC}=200$, seed $2026$.

---

## Four Parts (each looped over linkage $\in$ {Single, Complete, Average})

1. **SNR sweep** — $\mathrm{SNR}\in\{0.1,0.2,0.5,0.6,1.0,2.0,5.0\}$.
2. **$\rho$ sweep** — $\rho\in\{0.0,0.1,\dots,0.9\}$.
3. **$Q$ sweep** — $Q\in\{5,10,\dots,50\}$; $p_{ar}=M\times Q$ varies, $p_{\text{white}}=1000-p_{ar}$.
4. **$M$ sweep** — $M\in\{1,3,5,10,15,20,30\}$; $p_{ar}$, $p_{\text{white}}$, and $s=M$ all vary.

---

## Output Files

Written to `simulation_results/` (24 files = 12 scenario stems, one `.txt`+`.csv` pair each):

- `da_trex_mc_da_ar1_white_snr_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_white_rho_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_white_Q_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_white_M_{Single,Complete,Average}.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_da/demo_trex_da_05_mc_sim_bt_ar1_block_sweeps/demo_trex_da_05_mc_sim_bt_ar1_block_sweeps
```

---

## Real Results (excerpt: SNR sweep, Single linkage, TLARS)

| SNR | FDR | TPR |
|---|---|---|
| 0.10 | 0.0000 | 0.0160 |
| 0.50 | 0.0124 | 0.6850 |
| 1.00 | 0.0165 | 0.9040 |
| 2.00 | 0.0153 | 0.8880 |
| 5.00 | 0.0155 | 0.8970 |

---

## Interpretation

- FDR is **very tightly controlled** here (well below the already-lower $\mathrm{tFDR}=0.1$ target for this demo), and notably lower than Demo 04's block-only design at comparable SNR — with $975$ uncorrelated white-noise columns diluting the design, there is simply much less correlated "false-positive-prone" structure for the selector to latch onto compared to Demo 04's all-block design.
- TPR reaches a strong plateau around $0.88$–$0.90$ from $\mathrm{SNR}=1.0$ upward, similar in character to Demo 04 but slightly higher — suggesting the added white-noise columns do not meaningfully hurt detection of the true AR(1) blocks, consistent with T-Rex's general robustness to additional null variables.
- Compare directly against **Demo 04** (no white-noise columns, same block AR(1) core) to isolate the effect of a much larger overall $p$ with mostly-null structure versus a smaller, fully-block design.
- As in Demo 04, the $Q$/$M$/linkage sweeps are the key diagnostics for how block size, block count, and HAC linkage choice affect recovery in this white-noise-diluted setting.

---

**Last updated**: 2026-07-08
