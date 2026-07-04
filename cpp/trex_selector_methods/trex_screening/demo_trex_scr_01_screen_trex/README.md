# Demo 01: Screen-TRex (In-Memory)

## Purpose

Demonstrate basic Screen-TRex functionality: a single-run walkthrough of the **Ordinary** ($\Phi_j > 0.5$) and **Bootstrap-CI** screening thresholds, followed by a Monte Carlo SNR sweep comparing the two. In-memory variant — see Demo 02 for the memory-mapped counterpart.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

### Single-run demos (Demo 01 / 01b in the source)
- $n = 300$, $p = 1000$ (high-dimensional setting; the function also supports a low-dimensional $n=1000,p=300$ mode, not exercised by `main()`).
- Fixed true support $\{27, 149, 398, 4, 42\}$ ($s=5$), fixed coefficients $\beta_j = 1.0$.
- $\mathrm{SNR} = 5.0$.
- Demo 01 uses **Ordinary** screening; Demo 01b uses **Bootstrap-CI** screening, both on the same data.

### Monte Carlo SNR sweep (Demo 02 in the source)
- $n = 300$, $p = 1000$, random support of size $s = 10$ drawn fresh per trial.
- $\mathrm{SNR} \in \{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0\}$.
- $\mathrm{MC} = 500$ repetitions per (method, SNR) point.

---

## Control Parameters

```
K = 20                        # Dummy experiments per T-loop iteration
solver = TLARS
use_memory_mapping = false
R_boot = 1000                 # Bootstrap-CI replicates
ci_grid_step = 0.001
```

---

## Methods Compared

- **Screen-TRex Ordinary** ($\Phi_j > 0.5$ threshold)
- **Screen-TRex Bootstrap-CI**

---

## Output Files (expected)

Written to `simulation_results/` once run:

- `d01_screen_trex_mc_n300_p1000.txt` / `.csv` — MC SNR-sweep FDR/TPR/estimated-FDR table.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_scr_01_screen_trex
```

---

## Interpretation

- Both single-run demos print the selected indices, FDP/TPP, and estimated FDR for one dataset — useful to sanity-check that Screen-TRex runs end-to-end and returns a sensible selection before looking at aggregate MC behavior.
- In the MC sweep, expect both methods to show low TPR at very low SNR (0.01–0.1), rising toward high TPR as SNR increases, with FDR/estimated FDR reasonably close to each other — a large, persistent gap between realized FDR and the method's own estimated FDR would indicate the screening threshold's self-assessment is miscalibrated for this design.
- Bootstrap-CI is expected to be more conservative (and possibly slower, given $R_{\text{boot}}=1000$ replicates) than Ordinary; whether that translates into lower TPR at matched FDR is the main comparison to check once this demo has been run.

---

**Last updated**: 2026-07-04
