# Demo 04: DA-TRex (BT) on Block-Diagonal AR(1) Data

## Purpose

Introduce **`DAMethod::BT`** (binary-tree/dendrogram HAC clustering) on a clean block-diagonal AR(1) design, and study sensitivity to the HAC **linkage method** (Single/Complete/Average) across four sweep dimensions (SNR, $\rho$, block size $Q$, number of blocks $M$). Corresponds to R reference `demo_trex_da_05_bt_ar1_block_sweeps.R`.

---

## Data Generation Parameters (`dgp_ar1_block`)

- $p$ predictors partitioned into $M$ independent blocks of size $Q$ ($p = M \times Q$); within each block, an AR(1) recursion $X_{i,\text{col}} = \rho X_{i,\text{col}-1} + \sqrt{1-\rho^2}\eta$; zero correlation between blocks.
- **One active variable per block** (`OnePerBlock` support) — all $M$ blocks are active, so $s = M$.
- Base parameters: $n=150$, $M=5$, $Q=5$ ($p=25$, $s=5$), $\rho=0.7$, amplitude $1.0$, $\mathrm{tFDR}=0.2$, $K=20$, $\mathrm{MC}=200$, seed $2026$.

---

## Four Parts (each looped over linkage $\in$ {Single, Complete, Average})

1. **SNR sweep** — $\mathrm{SNR}\in\{0.1,0.2,0.5,0.6,1.0,2.0,5.0\}$.
2. **$\rho$ sweep** — $\rho\in\{0.0,0.1,\dots,0.9\}$.
3. **$Q$ sweep** — block size $Q\in\{5,10,\dots,50\}$; $p=M\times Q$ varies, $s=M=5$ fixed.
4. **$M$ sweep** — number of blocks $M\in\{1,3,5,10,15,20,30\}$; $p$ and $s=M$ both vary.

All 3 solvers (TLARS, TAFS, TOMP) run inside every linkage iteration.

---

## Output Files

Written to `simulation_results/` (24 files = 12 scenario stems, one `.txt`+`.csv` pair per {sweep} × {linkage} combination):

- `da_trex_mc_da_ar1_blk_snr_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blk_rho_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blk_Q_{Single,Complete,Average}.txt` / `.csv`
- `da_trex_mc_da_ar1_blk_M_{Single,Complete,Average}.txt` / `.csv`

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex_da/demo_trex_da_04_mc_sim_bt_ar1_block/demo_trex_da_04_mc_sim_bt_ar1_block
```

---

## Real Results (excerpt: SNR sweep, Single linkage, TLARS)

| SNR | FDR | TPR | Avg T |
|---|---|---|---|
| 0.10 | 0.0246 | 0.0260 | 11.20 |
| 0.50 | 0.1004 | 0.4430 | 22.06 |
| 1.00 | 0.0849 | 0.7640 | 27.88 |
| 2.00 | 0.0911 | 0.7650 | 26.55 |
| 5.00 | 0.0770 | 0.7930 | 28.14 |

---

## Interpretation

- FDR stays broadly near or somewhat above the $\mathrm{tFDR}=0.2$ target across the SNR grid, while TPR climbs from near-zero at $\mathrm{SNR}=0.1$ to a plateau around $0.76$–$0.79$ from $\mathrm{SNR}=1.0$ upward — the plateau (rather than reaching $\approx1.0$ as in the unblocked AR(1) case of Demo 01) suggests that recovering *every* one of the 5 active blocks consistently is harder than recovering isolated active variables, even at high SNR.
- The $Q$ and $M$ sweeps are the most novel diagnostics here: as block size $Q$ grows (more within-block noise variables per active signal) or the number of blocks $M$ grows (more simultaneous active blocks to find), expect the selection problem to become harder — check whether TPR degrades gracefully or drops sharply beyond some threshold.
- The linkage sweep answers a practical question for users of the `BT` method: does Single/Complete/Average linkage materially change FDR/TPR here, or is performance robust to that choice? Compare the three `_Single`/`_Complete`/`_Average` files for the same sweep to check.

---

**Last updated**: 2026-07-08
