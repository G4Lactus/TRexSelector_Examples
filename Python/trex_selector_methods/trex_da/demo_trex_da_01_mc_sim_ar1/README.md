# Demo 01: DA-TRex+AR1 on AR(1) Toeplitz Data (Python)

Monte-Carlo simulations for **DA-TRex+AR1** — the AR(1)-model dependency-aware
T-Rex selector — vs. the **classical (no-DA) T-Rex selector** on AR(1) Toeplitz
data (n=300, p=1000, s=10, amplitude 3.0, tFDR=0.2, K=20 random experiments,
MC=200 per grid point; solvers TLARS / TAFS / TOMP).

Python port of the C++ demo
[cpp/trex_selector_methods/trex_da/demo_trex_da_01_mc_sim_ar1](../../../../cpp/trex_selector_methods/trex_da/demo_trex_da_01_mc_sim_ar1/README.md)
— same DGP, control parameters, sweep axes, console summaries, and result-file
naming.

The greedy solvers use *exchangeable tie-breaking* (`exch_tie_alpha = 0.25` for
TAFS/TOMP, `0` for TLARS); TAFS additionally runs with its AFS correlation
parameter `rho_afs = 0.3` (`0` for TLARS/TOMP). With this fix, all three DA
solvers are FDR-controlled across the entire rho grid and both support
placements in the C++ reference results.

## Data-generating process

`dgp_ar1_snr` builds a Toeplitz design where `Sigma[j,k] = rho^|j-k|` via the
recursion `X[i,j] = rho * X[i,j-1] + sqrt(1 - rho^2) * eta[i,j]`, then
calibrates the noise to a target SNR
(`sigma^2 = Var(X @ beta) / SNR`). Supports are placed by the `CappedSpread`
policy (evenly spaced, gap capped at `max_gap=20`) or the `Random` policy
(uniform without replacement, redrawn per trial).

## What it runs

Parts are toggled by the `RUN_PART_*` flags at the top of the file (all `True`,
mirroring the cpp `if (true)` toggles).

- **Part 1**: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}` at rho=0.7,
  CappedSpread(max_gap=20) and Random support, each with the base (no-DA)
  T-Rex comparison rows (`TREX: <solver>`).
- **Part 2**: rho sweep `{0.0, 0.1, ..., 0.9}` at SNR=2.0, same two support
  runs with base comparison.
- **Part 3**: 2D gap x rho sweep; `gap_grid = {100, 50, 20, 15, 10, 5, 1}`,
  `rho_grid = {0.0, ..., 0.9}`, CappedSpread(gap) matrices plus a Random
  column, with the kappa boundary `kappa = ceil(log(0.02)/log(rho))` — the
  DA+AR1 correction-window half-width. When gap < kappa, active predictors
  fall inside each other's correction windows and TPP collapses.

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.py
```

The file inserts its own folder and the parent suite directory onto `sys.path`,
so the shared `dgp_generators` / `da_sim_common` modules import directly. Monte
Carlo trials run in parallel via `concurrent.futures.ProcessPoolExecutor`
(up to `min(6, os.cpu_count())` workers). Reduce `NUM_MC` to preview.

## Output files

Written to `simulation_results/` (`.txt` pretty-printed, `.csv` tidy long
format; the cpp suite splits into `data/` + `plots/`, the Python suite has no
plotting layer):

- `da_trex_mc_da_ar1_snr_capped.txt` / `.csv`
- `da_trex_mc_da_ar1_snr_random.txt` / `.csv`
- `da_trex_mc_da_ar1_rho_capped.txt` / `.csv`
- `da_trex_mc_da_ar1_rho_random.txt` / `.csv`
- `da_trex_mc_da_ar1_gap_rho.txt` / `.csv`

## What to expect

See the C++ counterpart README for the committed reference results: all three
DA solvers stay below tFDR=0.2 at every SNR while base T-Rex TLARS violates
from SNR=0.5 on (FDR ~0.29); at rho=0.9 the CappedSpread(max_gap=20) run is in
the collapse regime (gap 20 < kappa 38), while Random support keeps power with
controlled FDR. The gap x rho TPP heatmap shows the sharp kappa-boundary cliff
— a structural property of the AR1 correction window, not a solver artifact.

---

**Last updated**: 2026-07-21
