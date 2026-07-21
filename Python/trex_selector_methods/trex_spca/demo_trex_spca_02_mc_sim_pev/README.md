# Demo 02: Cumulative Explained Variance of T-Rex SPCA (Python)

## Purpose

While Demo 01 asks whether the loading support is FDR-controlled, this demo asks: **what proportion of
 variance does the selected support actually explain, and where does a non-sparse method's apparent advantage
 come from?**

It runs all four sweeps of **Fig. 3** of [[1]](#references) — cumulative explained variance under the
 paper's Definition 1 — and reports the conventional total-variance normalization alongside it so the two can
 be compared directly.
The pairing is the point: under one denominator ordinary PCA looks best, under the
 other it looks worst, and both statements are computed from the very same explained variance.

Python counterpart of
 [`cpp/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev/`](../../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md)
 (same DGP, sweeps, measurement definitions and output-file naming). Support indices here are **0-based**
 (Python convention).

The Definition-1 denominator is the **Signal + Mixed EV of the data** — the variance carried by the true
 active variables, fixed per dataset and shared by all methods. This convention reproduces the published
 panels' qualitative behaviour (curves *rising* with SNR from ≈40 %, ordinary PCA overshooting far past
 100 %, sparse methods saturating near 100 %). An earlier revision of the cpp demo instead normalized each
 method by its **own** Signal + Mixed EV ("Definition 1 as printed"), which is structurally ≥ 100 % for
 null-capturing methods and inverts the SNR trend — the published curves are unreachable under it.
 Digit-level coincidence with the paper's Fig. 3 remains impossible either way: the paper's exact
 normalization cannot be reconstructed from its text (see
 [the dedicated section below](#relation-to-the-papers-fig-3)).

---

## Data Generation Parameters (`generate_sparse_factor_model`)

The sparse $M$-factor model of [demo 01](../demo_trex_spca_01_mc_sim/README.md#data-generation-parameters-generate_sparse_factor_model),
 unchanged:

```math
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}
```

- $n = 50$, $p = 100$, $M = 3$ true factors, $\sigma_m \in \{5, 3, 1\}$.
- $p_1 = 10$ active loadings of value $0.9$ per factor, drawn without replacement from a **shared pool** of
   `overlap_pool_size` $= 30$ candidates, so factor supports partially coincide. Unlike demo 01 (which
   mirrors Fig. 2 at $p_1 = 5$), this demo mirrors Fig. 3 and therefore defaults to $p_1 = 10$.
- All methods see the same **center-only** $\boldsymbol{X}$ (covariance-PCA footing).

---

## Control Parameters

```text
tFDR = 0.10           # Target FDR (swept in Part 4)
lambda_2 = -1         # < 0 selects the ridge penalty by k-fold CV (0 = none, > 0 = fixed)
scaling = L2          # Per-component selector scaling mode
MC = 200              # Monte Carlo repetitions per grid point
base_seed = 42        # Per-grid-point seed band, derived from the swept value
```

As in demo 01, only the **data** is seeded deterministically; each trial's dummies come from hardware
 entropy (selector seed $-1$), which is what makes the Monte Carlo estimate valid. A re-run therefore
 reproduces results to within Monte Carlo noise, not exactly. One Python-specific detail: the cpp demo
 stores the per-grid-point seed offset in an `unsigned int`, which silently wraps the negative offsets the
 SNR sweep produces — the Python port mirrors that wrap explicitly (`% 2**32`), since numpy rejects
 negative seeds.

---

## Methods Compared

Identical to [demo 01](../demo_trex_spca_01_mc_sim/README.md#methods-compared): **OrdPCA**, **OraclePCA**,
 and the four T-Rex SPCA variants crossing `TENET`/`TENET_AUG` with `ActiveSet`/`Thresholded`.

The metrics are those of [demo 01](../demo_trex_spca_01_mc_sim/README.md#metrics) — classical EV, corrected
 to the adjusted EV of [[3]](#references), reported as `PEV` (total-variance share), `PEVsig`
 (Definition 1) and `PEVsigmix` (the method's own Sig + Mix part; OrdPCA's is the paper's
 "Ordinary PCA (Sig + Mix)" reference curve). As in demo 01, realized FDR is evaluated on PC1's support
 only, and FDR *control* is assumed for PC1 only: the overlapping factor supports make later components'
 ground truth ambiguous.

---

## The Sweeps

Each part varies **one** axis and holds the rest at their defaults (#PCs $= 3$, SNR $= 0$ dB, $p_1 = 10$,
 $\mathrm{tFDR} = 0.10$):

| Part | Axis | Grid |
| --- | --- | --- |
| 1 | Number of **extracted** PCs | $\{1, 3, 5, 10, 20, 30, 40, 49\}$ |
| 2 | SNR (dB) | $\{-10, -7, -5, -3, 0, 3, 5, 7, 10\}$ |
| 3 | Number of true active loadings $p_1$ | $\{1, 5, 10, 15, 20, 25, 30\}$ |
| 4 | Target FDR | $\{0.01, 0.05, 0.10, 0.15, \ldots, 0.50\}$ |
| 5 | Realized FDR / TPR **heatmaps** over the full target-FDR $\times$ SNR grid | $11 \times 9$ points, PC1 support, two selector responses |

Notes on the grids:

- **Part 1 varies only how many components are *extracted*.** The data keeps its $M = 3$ true factors
   throughout, so every component beyond the third is null by construction and can contribute nothing but
   Null EV. That is exactly what makes the panel diagnostic.
- **Part 1 stops at 49, not at 50.** $\boldsymbol{X}$ is centered, so $\operatorname{rank}(\boldsymbol{X})
   \leq n - 1 = 49$ and a 50th component carries no variance at all — asking the per-PC selector to regress
   on that null direction throws outright, because its response vector is constant. 49 is the largest
   attainable PC count on an $n = 50$ design, not a safety margin.
- **Part 3 is bounded above by the overlap pool.** Each factor draws its support from a shared pool of 30
   candidates, so $p_1 \leq 30$.
- **Only the T-Rex variants respond to Part 4.** The two PCA baselines never see `tFDR`; they are reported
   as flat reference lines.
- **Parts 2–4 are displayed as cumulative PEV only.** Their extraction is held at the $M = 3$ true factors,
   so what they show is the cumulative PEV of PCs 1–3 under both normalizations — the tables carry nothing
   else. Support-recovery metrics belong to demo 01's sweep; here they remain in the CSVs as the data
   record but are not displayed.
- **Part 5 is the support-recovery counterpart of Part 4**, run over the full 2D grid instead of the 0 dB
   slice, and only for PC1. Each trial runs the per-PC T-Rex+GVS(EN) selector
   (`tsn.TRexGVSSelector`, the same configuration `tsn.TRexSPCASelector` uses internally) twice on the same
   data: once on the **plug-in response** $y = \boldsymbol{X} v_1$ (what T-Rex SPCA actually regresses on)
   and once on the **oracle response** $y = z_1$, the true factor scores. The pairing isolates the mechanism
   behind the FDR overshoot at loose targets — see the caveat section below.

---

## Output Files

Written to `simulation_results/data/`, one pair per sweep:

- `demo_trex_spca_02_mc_sim_pev_pcs.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_snr.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_loadings.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_tfdr.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_fdr_heatmap.csv` — Part 5 (`response,metric,tfdr,snr_db,value`)

The `pcs` table carries FDR, TPR, `PEV` (total-variance normalization), `PEVsig` (Definition 1) and
 `PEVsigmix` (the method's Sig + Mix part) per method; the `snr`, `loadings` and `tfdr` **tables are
 restricted to the cumulative-PEV rows** (PCs 1–3). Every CSV still records all five metrics.

Figures (PNG + PDF) go to `simulation_results/plots/`, produced by `./generate_plots.sh` via the shared
 suite plotting module [`../trex_spca_plt_utils.py`](../trex_spca_plt_utils.py), which combines all four
 sweeps into one panel figure per normalization:

- `demo_trex_spca_02_mc_sim_pev_signal.png` — cumulative PEV under Definition 1.
- `demo_trex_spca_02_mc_sim_pev_total.png` — cumulative PEV as a share of the total variance.
- `demo_trex_spca_02_mc_sim_pev_fdr_heatmap.png` — Part 5: realized FDR / TPR heatmaps, both responses.

> **No committed results yet**: this Python port has been verified end-to-end with pilot runs (reduced MC
> counts and grids), but the full 200-trial sweeps have not been run here — `simulation_results/` is
> intentionally absent until a full run is committed. For the full-quality numbers and the result
> discussion, see the
> [cpp demo-02 README](../../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md);
> the pilot runs reproduce its qualitative behaviour (see the summary below).

---

## Running the Demo

```bash
cd Python/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev
python demo_trex_spca_02_mc_sim_pev.py            # default min(6, cores) workers
python demo_trex_spca_02_mc_sim_pev.py 8          # 8 worker cores
./generate_plots.sh                               # render the figures from the saved CSVs
```

`TREX_SPCA_NUM_MC=<k>` overrides the Monte Carlo count for quick pilot runs (e.g. `TREX_SPCA_NUM_MC=10`),
and `TREX_SPCA_PARTS` selects which parts run (default `12345`; e.g. `TREX_SPCA_PARTS=5` reruns only the
FDR/TPR heatmap sweep) — the same environment switches as the cpp executable.

Run the demo inside an environment that carries `trex_selector_neo` (e.g. the `trex-dev` conda environment
or the repo `.venv`); it parallelizes trials via `concurrent.futures.ProcessPoolExecutor` (spawn), so run it
**as a real `.py` file**, not via `python -` / piped stdin. The optional first CLI argument sets the number
of worker cores.

The PC-count sweep dominates the runtime: extracting 49 components means 49 independent T-Rex selections
per trial, so Part 1 alone costs several times what the other parts cost together — expect a full run to
take substantially longer than in C++ (Python drives the same C++ core, but per-trial process dispatch and
the per-PC CV add overhead).

---

## What the cpp reference finds

The full-quality cpp run (200 MC per grid point) establishes the reading this port reproduces
 qualitatively in pilot runs:

- **PC-count sweep (panel a)**: ordinary PCA climbs past the 100 % Definition-1 line at about five PCs and
   reaches ≈165 % at 49 — reporting two thirds more variance than the true active variables carry — while
   the FDR-controlled variants flatten out by three-to-five PCs and stay flat, declining to claim variance
   that is not theirs. OraclePCA overshoots too: knowing the support *size* does not protect against
   thresholding noise components.
- **SNR sweep (panel b)** rises with SNR as in the published Fig. 3(b); ordinary PCA's head start at low
   SNR is exactly its Null EV, as its dashed Sig + Mix decomposition shows.
- **Loadings sweep (panel c)** is non-monotone with a minimum around $p_1 = 15$–$20$ — the middle of the
   range is where support selection actually matters.
- **tFDR sweep (panel d)**: T-Rex pays visibly only at the strictest target ($\mathrm{tFDR} = 0.01$), then
   climbs onto a gentle plateau — loosening the target buys captured variance slowly.
- **Heatmaps (part 5)**: with the plug-in response, realized FDR holds up to the targets this suite
   actually uses ($\leq 0.10$) and breaks systematically from $\mathrm{tFDR} = 0.15$ upward; with the
   oracle response the overshoot largely vanishes — confirming it is a property of scoring a dense-truth
   regression against a sparse truth, not of the selector. **T-Rex SPCA controls the FDR of the regression
   it is handed, and that transfers to the factor model only insofar as the plug-in PC is a faithful proxy
   for the factor.**

### Relation to the paper's Fig. 3

The paper's text does not pin down the normalization behind its plotted curves — no constructible
 denominator reproduces all published anchors, so digit-level agreement with Fig. 3 is not attainable from
 the paper alone. (A self-normalized reading of Definition 1 is ruled out entirely: it is $\geq 100\%$ by
 construction for null-capturing methods and inverts the SNR trend.) This demo therefore uses the fixed
 active-column denominator throughout, which reproduces the figure's qualitative content: curves rising
 with SNR, ordinary PCA's Null-EV overshoot past $100\%$, the sparse methods saturating near $100\%$, and
 the ranking inversion between the two normalizations.

---

## References

1. Machkour, J., Breloy, A., Muma, M., Palomar, D. P., & Pascal, F., "Sparse PCA with False Discovery Rate Controlled
   Variable Selection.", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2024,
   pp. 9716–9720, IEEE.
   [DOI-Link](https://doi.org/10.1109/ICASSP48485.2024.10448237)
2. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Grouped Variable Selection
   in High-Dimensional Linear Models using the T-Knock Filter.", European Signal Processing Conference (EUSIPCO), 2022,
    pp. 892–896, EURASIP.
    [DOI-Link](https://doi.org/10.23919/EUSIPCO55093.2022.9909883)
3. Zou, H., Hastie, T., & Tibshirani, R., "Sparse Principal Component Analysis.", Journal of Computational
   and Graphical Statistics, vol. 15, no. 2, 2006, pp. 265–286, Taylor & Francis.
   [DOI-Link](https://doi.org/10.1198/106186006X113430)

---

**Last updated**: 2026-07-21
