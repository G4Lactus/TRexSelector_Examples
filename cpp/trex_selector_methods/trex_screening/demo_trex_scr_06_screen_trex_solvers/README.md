# Demo 06: Screen-TRex Solver Comparison

## Purpose

Benchmark whether the underlying T-Rex solver backend — **TLARS**, **TAFS** ($\rho_{\text{afs}}=0.3$), or **TOMP** — changes Screen-TRex's FDR/TPR behavior, under both the Ordinary and Bootstrap-CI screening thresholds.

> **Status**: `simulation_results/` is currently empty — this demo has not yet been run in this checkout.

---

## Data Generation Parameters

- $n=300$, $p=1000$ (high-dimensional setting; the underlying functions also support a low-dimensional $n=1000,p=300$ mode, not exercised by `main()`).
- Random support of size $s=10$ drawn fresh per trial.
- $\mathrm{SNR} \in \{0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0\}$.
- $\mathrm{MC} = 500$ repetitions per (solver/method, SNR) point.

---

## Control Parameters

```
K = 20
tloop_stagnation_stop = true
tloop_max_stagnant_steps = 5
solver-specific:
  TLARS:     (no extra parameters)
  TAFS-0.3:  rho_afs = 0.3
  TOMP:      (no extra parameters)
```

---

## Two Parts

1. **Solver comparison, Ordinary only** — 3 series: `TLARS (Ordinary)`, `TAFS-0.3 (Ordinary)`, `TOMP (Ordinary)`.
2. **Solver × method comparison** — 6 series: the same 3 solvers, each under both `(Ordinary)` and `(Bootstrap)` Screen-TRex.

---

## Output Files (expected)

Written to `simulation_results/` once run:

- `d06_screen_trex_solvers_mc_snr_n300_p1000_s10.txt` / `.csv` — 3-series solver comparison (Ordinary only).
- `d06_screen_trex_solver_method_mc_snr_n300_p1000_s10.txt` / `.csv` — 6-series solver × method comparison.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_scr_06_screen_trex_solvers
```

---

## Interpretation

- If Screen-TRex's screening decision depends mainly on the dummy-based voting proportion $\Phi_j$ rather than on solver-specific path details, expect **TLARS**, **TAFS-0.3**, and **TOMP** to produce broadly similar FDR/TPR curves — with differences concentrated in runtime rather than accuracy, mirroring the solver-equivalence pattern seen in the classical T-Rex demos (see `../../trex/demo_trex_02_mc_sim_fixed_support/README.md`).
- The 6-series part is the more informative comparison: check whether the **Ordinary vs. Bootstrap-CI** gap (if any) is consistent across all three solvers, or whether one solver interacts differently with the bootstrap thresholding step.
- `TAFS-0.3`'s extra regularization parameter ($\rho_{\text{afs}}=0.3$) is the main place a solver-specific difference could plausibly show up, since it changes the T-loop's adaptive forward-stagewise behavior rather than just the path algorithm.

---

**Last updated**: 2026-07-04
