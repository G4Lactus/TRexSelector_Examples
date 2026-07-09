# Demo 01: Screen-TRex (In-Memory)

## Purpose

Basic in-memory Screen-TRex variable screening: two single-run walkthroughs of the **Ordinary** (`Phi_j > 0.5`) and **Bootstrap-CI** thresholding rules on the same fixed-support dataset, followed by a Monte Carlo SNR sweep comparing the two on FDR, TPR, and the procedure's self-reported Estimated FDR. In-memory variant — see Demo 02 for the memory-mapped counterpart. Mirrors `cpp/trex_selector_methods/trex_screening/demo_trex_scr_01_screen_trex/`.

## Data-generating process

`dgp_iid_snr` — i.i.d. standard-Gaussian design with SNR-calibrated noise, `sigma = sqrt(||X beta||^2 / (n * SNR))`.

- **Single runs (Parts 1–2)**: `n = 300`, `p = 1000` (high-dimensional, `p > n`), fixed 0-based true support `[27, 149, 398, 4, 42]` (`s = 5`), coefficients `beta_j = 1`, `SNR = 5.0`, data seed 58, selector seed 42.
- **MC sweep (Part 3)**: `n = 300`, `p = 1000`, random support of size `s = 10` drawn fresh per trial, `SNR in {0.01, 0.1, 0.2, 0.5, 0.6, 1.0, 2.0, 5.0}`, `NUM_MC = 60` (the C++ demo uses 500; reduced for a ~1-minute runtime).

## Methods compared

- **Screen-TRex Ordinary** — `make_screen_ctrl(trex_method="TREX", bootstrap=False)`
- **Screen-TRex Bootstrap-CI** — same, with `bootstrap=True` (library defaults for the bootstrap replicates / CI grid step)

Shared solver controls via `make_scr_trex_ctrl()`: `solver = "TLARS"`, `K = 20`, `use_memory_mapping = False`. Selection is driven by `ts.TRexScreeningSelector(X, y, screen_control=..., trex_control=..., seed=42)`; `.select()` returns a `ScreenTRexSelectionResult` exposing `estimated_FDR`, `estimated_correlation`, `used_bootstrap`, and `confidence_level`, with the chosen variables in `sel.selected_indices` (0-based NumPy array).

## Parts

- **Part 1** — Single-run Ordinary Screen-TRex: selected indices (0-based), FDP/TPP against the known support, Estimated FDR.
- **Part 2** — Single-run Bootstrap-CI Screen-TRex on the same data.
- **Part 3** — Monte Carlo SNR sweep, Ordinary vs. Bootstrap-CI, averaging FDR / TPR / Estimated FDR over 60 runs per (method, SNR) point.

## Running

```bash
python Python/trex_selector_methods/trex_screening/demo_trex_scr_01_screen/demo_trex_scr_01_screen.py
```

A small `sys.path` bootstrap at the top of the script inserts the demo directory and its parent (suite) directory, so the shared `trex_scr_common.py` module resolves regardless of the working directory. The `RUN_PART_*` flags at the top toggle individual parts. Parts 1–2 print to the console only; Part 3 additionally writes `simulation_results/d01_screen_trex_mc_n300_p1000.txt` (summary table) and `.csv` (tidy long format). Until the demo is run, `simulation_results/` holds only a `.gitkeep`. Results are statistically comparable to the C++ reference, not bit-identical (different RNG streams).

## Interpretation

- The single runs sanity-check that Screen-TRex runs end-to-end and returns a sensible selection at a comfortable SNR of 5 before looking at aggregate MC behavior.
- Expect Ordinary FDR to fall from its high-noise value at SNR 0.01 toward a small residual for SNR >= 1, with TPR rising from near zero to 1.00 at SNR 5; the Estimated FDR should track the realized FDR reasonably, typically sitting somewhat above it at high SNR.
- Bootstrap-CI is the more conservative rule throughout: lower FDR but a clear TPR cost, and it is slower per run due to the bootstrap replicates.

**Last updated**: 2026-07-08
