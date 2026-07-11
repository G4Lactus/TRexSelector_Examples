# T-Rex+GVS: Grouped Variable Selection — Python Demos

Python ports of the **Grouped Variable Selection T-Rex (T-Rex+GVS)** demos, for
designs where the active variables cluster into groups. They mirror the R suite
in `R/trex_selector_methods/trex_gvs/` and the C++ suite in
`cpp/trex_selector_methods/trex_gvs/`, using the `trex_selector_neo`
`TRexGVSSelector`.

The folder layout mirrors the C++/R suites one-to-one: one subfolder per demo
(each with its own `README.md` and `simulation_results/`), plus two suite-level
shared modules.

```txt
trex_gvs/
  ├── README.md
  ├── trex_gvs_dgps.py        ┐  suite-level shared modules
  ├── gvs_sim_common.py       ┘  (imported by every demo)
  ├── demo_trex_gvs_01/
  │   ├── demo_trex_gvs_01.py
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_gvs_02/
  │   └── ...
  └── demo_trex_gvs_08/
```

Shared infrastructure:

- [trex_gvs_dgps.py](trex_gvs_dgps.py) — all scenario data-generating processes
  in one module (mirroring the C++ `trex_gvs_dgps.hpp`): `dgp_hastie_snr`,
  `dgp_scattered_grouped_snr`, `dgp_mixed_blocks_snr`, `dgp_neg_traps_snr`,
  `dgp_t3_equi_blocks_snr`, `dgp_ar1_blocks_snr`, `dgp_arma_blocks_snr`
  (+ the `_sim_arma` helper), and `dgp_block_equicorr`. All are SNR-calibrated
  (`sigma_y = sqrt(var(signal) / snr)`, `ddof=1`; the block benchmark uses the
  biased `ddof=0` denominator, matching the cpp DGP). Support / group indices
  are 0-based.
- [gvs_sim_common.py](gvs_sim_common.py) — the `build_gvs_selector` adapter
  (EN / EN+AUG / IEN, HAC or oracle `prior_groups`), a generic parallel Monte
  Carlo runner `run_mc_gvs`, and the `print_table` / `print_matrix` /
  `print_single_run_result` printers.

Each demo begins with a **Part 1 single run** (data/structure diagnostics plus
one EN / EN+AUG / IEN run) and then offers Monte Carlo sweeps that compare the
three dummy-generation methods **EN (TENET) / EN+AUG (TENET_AUG) / IEN**,
sweeping SNR, a correlation knob, and a 2-D SNR × correlation grid. Typical
fixed settings: `n = 200`, `p = 500`, `tFDR = 0.1`, `K = 20`, 200 MC runs
(demo 08 uses 500).

Dependencies: `trex_selector_neo`, `numpy`. No plotting dependency — the Python
Part 1 prints structural diagnostics rather than rendering the R suite's Plotly
correlation heatmap.

---

## Demos

| # | Scenario | Sweeps | DGP (in `trex_gvs_dgps.py`) | cpp counterpart |
|---|---|---|---|---|
| [01](demo_trex_gvs_01/) | Hastie: 3 equicorrelated active groups of 50 (all active, s = 150) | SNR, rho, 2-D, **lambda2_method**, **hc_linkage** | `dgp_hastie_snr` | `demo_trex_gvs_01_mc_sim_hastie_en_blocks` |
| [02](demo_trex_gvs_02/) | Scattered-grouped: same 3 groups, members permuted across columns | SNR, rho, 2-D | `dgp_scattered_grouped_snr` | `demo_trex_gvs_02_mc_sim_scattered_blocks` |
| [03](demo_trex_gvs_03/) | Mixed blocks: 3 active blocks (20/50/80) + 1 trap (65), random order, noise gaps | SNR, rho, 2-D — **two presets** (default sd_x; fixed rho = 0.75) | `dgp_mixed_blocks_snr` | `demo_trex_gvs_03_mc_sim_mixed_blocks` |
| [04](demo_trex_gvs_04/) | Negative traps: active +Z/−Z group + two inactive +Z/−Z traps (s = 100) | SNR, rho, 2-D | `dgp_neg_traps_snr` | `demo_trex_gvs_04_mc_sim_neg_traps` |
| [05](demo_trex_gvs_05/) | Heavy-tailed blocks: demo-03 geometry with Student-t(3) factors/noise | SNR (fixed rho = 0.75), rho, 2-D | `dgp_t3_equi_blocks_snr` | `demo_trex_gvs_05_mc_sim_t3_blocks` |
| [06](demo_trex_gvs_06/) | AR(1) blocks: demo-03 geometry with within-block Cor = rho^\|j−k\| | SNR (fixed rho = 0.85), rho, 2-D | `dgp_ar1_blocks_snr` | `demo_trex_gvs_06_mc_sim_ar1_blocks` |
| [07](demo_trex_gvs_07/) | ARMA blocks: AR(1) / MA(3) / ARMA(2,1) active blocks + AR(1) trap | SNR (fixed ar_coef = 0.8), **ar_coef**, 2-D | `dgp_arma_blocks_snr` | `demo_trex_gvs_07_mc_sim_arma_blocks` |
| [08](demo_trex_gvs_08/) | Block benchmark: 100 blocks × 5, 10 active; HAC-discovered vs oracle groups | 4 methods (M1–M4) × 3 truth scenarios over a (rho × snr) grid | `dgp_block_equicorr` | `demo_trex_gvs_08_mc_sim_block_bench` |

Each demo folder has its own `README.md`. Every demo exposes `RUN_PART_*` flags
at the top: **Part 1 (single run) runs by default; the heavy Monte Carlo parts
are opt-in** (set the flags to `True`). This differs from the R suite, where all
parts default on. Output is console-only.

### Demo 08 notes (block benchmark)

Demo 08 compares four variants — **M1** = EN + HAC-discovered groups, **M2** =
EN + oracle groups, **M3** = IEN + HAC, **M4** = IEN + oracle — across three
within-block truth scenarios (**INDIVIDUAL**: 1 active var/block;
**REPRESENTATIVE**: 2–3/block; **WHOLE**: all 5/block). Oracle groups are
supplied via the GVS control's `prior_groups` field (0-based contiguous labels);
HAC variants discover groups from the data (`corr_max = 0.5`, single linkage). It
reports coordinate TPR/FDR, block-level FDR, and a scenario-specific block metric
(`blk_hit` / `full_blk`), plus `T_stop`, `purity`, and `M_found`.

**Grouping diagnostics.** The `GVSSelectionResult` returned by `select()` exposes
`max_clusters` (→ `M_found`, the number of clusters formed) and `groups_vec`
(per-variable cluster labels), so demo 08 computes the real block-`purity` rate
(per the cpp benchmark) and `M_found` for the HAC methods as well as the oracle
ones. This surfaces more than the R demo, whose R6 API exposes only
`$max_clusters` and so reports `purity` as `1.0`/`-` rather than the genuine
per-block rate.

---

## Running

```bash
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_01/demo_trex_gvs_01.py       # Part 1 only (default)
.venv/bin/python Python/trex_selector_methods/trex_gvs/demo_trex_gvs_08/demo_trex_gvs_08.py 8     # 8 worker cores
```

Each demo inserts both its own directory and the parent suite directory onto
`sys.path`, so the shared `trex_gvs_dgps` / `gvs_sim_common` modules resolve from
the nested location and are re-importable by the `spawn`-launched Monte Carlo
workers — so demos must be run **as files**, not piped via `python -`. The
optional first CLI argument sets the number of parallel worker cores (default
`min(6, os.cpu_count())`). Enable the heavy sweeps by setting the `RUN_PART_*`
flags to `True`; reduce `num_MC` / grid sizes for a quick smoke run.

Indices are **0-based** (Python convention); the R counterparts are 1-based.

---

## Recorded results

Manually-recorded result summaries live in each demo's `simulation_results/`
folder, e.g.
[demo_trex_gvs_01/simulation_results/data/demo_trex_gvs_01_hastie_models.md](demo_trex_gvs_01/simulation_results/data/demo_trex_gvs_01_hastie_models.md).

---

**Last updated**: 2026-07-08
