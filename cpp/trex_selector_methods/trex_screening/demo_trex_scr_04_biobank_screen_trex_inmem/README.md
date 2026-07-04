# Demo 04: Biobank Screen-TRex (In-Memory)

## Purpose

Demonstrate the **Biobank Screen-TRex** workflow ("Algorithm 1"): screening **multiple phenotypes against one shared design matrix** $\boldsymbol{X}$, with adaptive per-phenotype routing between Screen-TRex Ordinary, Screen-TRex Bootstrap-CI, and a classical T-Rex fallback. In-memory variant — see Demo 05 for the memory-mapped counterpart.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

### Single phenotype
- $n=300$, $p=1000$, true support $\{4, 27, 42, 149, 398\}$ ($s=5$), $\beta_j = 1.0$, $\mathrm{SNR}=1.0$.

### Multiple phenotypes ($q=5$, shared $\boldsymbol{X}$)
- $n=300$, $p=1000$, fixed `X_seed` so all 5 phenotypes share the exact same design matrix.
- Per-phenotype support sizes $\{5, 6, 5, 7, 5\}$, each with $\beta_j = 1.0$ on its own random support.
- Per-phenotype $\mathrm{SNR} \in \{1.0, 2.0, 5.0, 10.0, 20.0\}$ (one SNR level per phenotype).

### Monte Carlo — single phenotype
- $n=300$, $p=1000$, random support $s=10$ per trial, $\mathrm{SNR} \in \{0.1, 0.5, 1.0, 2.0, 5.0\}$, $\mathrm{MC}=500$.

### Monte Carlo — multiple phenotypes
- $n=300$, $p=1000$, $q=5$ phenotypes per trial, support $s=5$ per phenotype, same SNR grid, $\mathrm{MC}=500$ (so $q \times \mathrm{MC} = 2500$ total phenotype screenings per SNR point).

---

## Control Parameters

```
lower_bound_FDR = 0.05
upper_bound_FDR = 0.15
target_FDR_trex = 0.10        # classical T-Rex fallback target
R_boot = 1000
ci_grid_step = 0.001
K = 20                        # T-Rex fallback dummy experiments
solver = TLARS
use_memory_mapping = false
```

---

## Methods (Adaptive Routing)

1. **Screen-TRex (Ordinary)** — tried first.
2. **Screen-TRex (Bootstrap-CI)** — tried if Ordinary's estimated FDR falls outside $[0.05, 0.15]$.
3. **T-Rex (fallback)** — classical T-Rex at fixed target FDR $0.10$, used if Bootstrap-CI also fails.

---

## Output Files (expected)

Written to `simulation_results/` once run:

- `d04_biobank_screen_trex_mc_snr_n300_p1000_s10_inmem.txt` / `.csv` — single-phenotype MC SNR sweep.
- `d04_biobank_screen_trex_mc_multi_n300_p1000_q5_s5_inmem.txt` / `.csv` — multi-phenotype MC SNR sweep.

Both include a **Usage %** column per method per SNR level. The single-phenotype and multi-phenotype demo functions (without "MC" in their name) print single-run/table diagnostics to the console only.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_scr_04_biobank_screen_trex_inmem
```

---

## Interpretation

- The single-phenotype and multi-phenotype (non-MC) demos are best used to sanity-check the routing logic on one dataset: check `method_used` per phenotype and confirm it is consistent with the reported `estimated_FDR` falling inside or outside the $[0.05, 0.15]$ window.
- In the MC sweeps, expect **Usage %** for **Ordinary** to increase with SNR (easier signal needs no fallback), with **Bootstrap-CI** and the **T-Rex fallback** picking up more of the routing share at lower SNR where Ordinary's estimated FDR is more likely to drift outside the acceptable window.
- Conditional FDP/TPP (computed only over trials routed to a given method) should be compared with the **unconditional** estimated FDR — since these are conditioned on different, possibly non-random subsets of trials, direct comparison across methods should account for the Usage % context rather than reading FDP/TPP columns in isolation.
- Comparing the single- vs. multi-phenotype MC results is informative: since all phenotypes in the multi-phenotype scenario share the same $\boldsymbol{X}$, any systematic difference from the single-phenotype case would point to cross-phenotype interaction effects (e.g. shared dummy/D-matrix reuse) rather than genuine per-phenotype statistical differences.

---

**Last updated**: 2026-07-04
