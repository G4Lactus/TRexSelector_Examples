# Demo 03: Monte Carlo Simulation with Variable Support Indices

## Purpose

Investigate T-Rex selector performance (FDR, TPR, average dummy multiplier $L$, average stopping time $T$) when the **positions** of the active variables are redrawn on every Monte Carlo trial, while the **number** of active variables stays fixed. This isolates the effect of support location (rather than support size) and complements the fixed-support study in Demo 02.

---

## Data Generation Parameters

- **Sample size**: $n = 300$
- **Number of features**: $p = 1000$
- **True support cardinality**: fixed at $s = 10$ on every trial. The support **size does not vary**; only the **indices** change — each trial shuffles $\{0, \ldots, p-1\}$ with `std::mt19937(seed + 500000)`, keeps the first 10, and sorts them.
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = false`)
- **SNR grid**: a **21-point grid** — $\{0.1, 0.2, \ldots, 2.0\}$ (step $0.1$) plus $5.0$
- **Monte Carlo repetitions**: `num_MC = 10` trials per solver × SNR level
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, \sigma^2 I_n)$, Normal predictors and Normal noise

`main()` runs only the high-dimensional configuration ($n = 300$, $p = 1000$).

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max dummies L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
dummy_distribution = Normal      # Dummy predictors drawn from N(0,1)
lloop_strategy = HCONCAT         # Horizontally concatenated dummy columns
tloop_stagnation_stop = true     # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 5     # Stagnation window (5 consecutive unchanged iterations)
tFDR = 0.1                       # Target FDR control level
```

The MC loop is parallelized with OpenMP (`omp_set_num_threads(6)`).

---

## Solvers Compared

The same **14** T-Rex base solvers as Demo 02, via the shared `make_default_solvers_to_test()` list (TLARS, TLASSO, TENET, TSTAGEWISE, TSTEPWISE, TOMP, TGP, TACGP, TMP, TAFS_rho_0.3, TAFS_rho_1.0, TNCGMP_v1, TNCGMP_v0, TOOLS). Unlike Demo 02, this demo also records the averaged dummy multiplier $L$ and stopping time $T$ per solver × SNR.

---

## Output Files

Both files are written to `simulation_results/data/`. The stem uses `tloop_max_stagnant_steps = 5` (hence `stagnation_window_5`), prefixed with `demo_trex_03_mc_sim_variable_support_`.

### Main Result File
**`demo_trex_03_mc_sim_variable_support_trex_results_n300_p1000_stagnation_window_5.txt`**

Aligned table with four metric rows per solver — FDR, TPR, Avg L, Avg T — across the 21 SNR columns.

### Tidy-Format CSV
**`demo_trex_03_mc_sim_variable_support_trex_results_n300_p1000_stagnation_window_5.csv`**

Long/stacked format, header column order **`solver,metric,snr,value`**, with `FDR`, `TPR`, `AvgL`, and `AvgT` rows per solver and SNR:
```
solver,metric,snr,value
TLARS,FDR,0.100000,...
TLARS,TPR,0.100000,...
TLARS,AvgL,0.100000,...
TLARS,AvgT,0.100000,...
...
```

---

## Running the Demo

```bash
./build/debug/bin/trex_selector_methods/trex/demo_trex_03_mc_sim_variable_support/demo_trex_03_mc_sim_variable_support
```

### Note on `TNCGMP_v0` memory messages (benign)

With older library builds, the `TNCGMP_v0` (LineSearch / Matching Pursuit) block prints lines such as:

```text
[TNCGMP_LineSearch] [WARNING] LineSearch worst-case beta path is ~1.34 GiB (3000 x 60000 doubles); rely on T_stop early stopping.
```

This is **not** an error and does **not** indicate a failed run — the reported results are correct. It is a *worst-case* memory estimate: the LineSearch variant allows atom re-selection and therefore uses a uniquely large step ceiling (`20 · p_tot`, vs `8 · min(...)` for the OMP-family solvers). The beta-path buffer is allocated **incrementally** and bounded by `T_stop` / stagnation stopping, so the full worst case is virtually never reached. Only `v0` reports it; `TNCGMP_v1` (FullyCorrective) does not. The message also fires once per solver construction (i.e. per MC trial), which is why it can repeat.

Upstream this has been reclassified from an unconditional `logWarning` to a verbose-only `logInfo`, so recent library builds are silent here unless verbosity is enabled.

---

## Key Questions Addressed

1. **Is FDR control robust to the location of the active set?**
   - Expected: FDR remains $\leq$ tFDR regardless of which 10 indices are active.

2. **How do the averaged $L$ and $T$ behave across SNR?**
   - The Avg L and Avg T rows expose how the dummy multiplier and stopping time adapt as signal strength rises.

3. **How does TPR climb across the fine SNR grid?**
   - The 21-point grid gives a smooth power curve from weak to strong signal.

---

## Interpretation Guide

**What to look for:**
- **FDR stability**: Should remain $\leq$ tFDR even as support indices change per trial.
- **TPR curve**: Smooth increase across the fine SNR grid.
- **Avg L / Avg T**: Diagnostic behaviour of the L-loop and T-loop across SNR.

**Comparison with Demo 02:**
- Demo 02 (fixed support indices) is the baseline.
- Demo 03 (support indices redrawn per trial, fixed cardinality) tests robustness to support location. Note the smaller MC count (`num_MC = 10`) — averages are noisier than Demo 02's 200 trials.

---

**Last updated**: 2026-07-08
