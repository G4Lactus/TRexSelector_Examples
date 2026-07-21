# T-Rex SPCA: Sparse Principal Component Analysis — Python Demos

## Overview

This folder introduces the sparse principal component version of the T-Rex selector, the **T-Rex SPCA**
 according to [[1]](#references).
A sparse principal component analysis method built on the T-Rex+GVS elastic-net machinery
 (see [../trex_gvs/](../trex_gvs/README.md)) according to [[2]](#references).

Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each
 principal component: the T-Rex+GVS selector is applied per component to control the false discovery rate of
 the estimated loading support. The result is a PCA whose components are supported on a small,
 FDR-controlled set of variables rather than on all $p$ of them.

The demos are the Python counterparts of the C++ suite in
 [`cpp/trex_selector_methods/trex_spca/`](../../../cpp/trex_selector_methods/trex_spca/README.md)
 (and the R suite in `R/trex_selector_methods/trex_spca/`), using the `trex_selector_neo` root-level
 `tsn.TRexSPCASelector`. They are designed to help users understand:

1. whether the FDR of the estimated loading support is actually held at the target,
2. what that control costs in explained variance, measured against ordinary PCA and an oracle that is told
    the true support size,
3. how the elastic-net solver (`TENET` vs. `TENET_AUG`) and the loading-assembly mode (`ActiveSet` vs.
    `Thresholded`) affect the outcome,
4. and why the *choice of denominator* for the explained variance decides which method looks best.

---

## What this folder covers

- A **sparse $M$-factor data-generating process** where each column of the loading matrix has exactly $p_1$
   nonzero entries drawn from a shared overlap pool, so factor supports can partially coincide.
- Two **T-Rex SPCA** solver variants crossed with two loading-assembly modes, benchmarked against **ordinary
   PCA** (no sparsity) and **oracle-thresholded PCA** (knows the true support size).
- A sweep over the signal-to-noise ratio **in decibels**, the scale this suite reports on throughout
   (demo 01), plus sweeps over the extracted PC count, the number of true active loadings and the target FDR
   (demo 02).

---

## Statistical model and targets

### The factor model

```math
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E},
\qquad
\boldsymbol{Z}_{\cdot,m} \sim \mathcal{N}(\boldsymbol{0}, \sigma_m^2 \boldsymbol{I}_n).
```

### Notation

- $n$: number of observations, $p$: number of variables, $M$: number of latent factors.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$: centered observation matrix.
- $\boldsymbol{V} \in \mathbb{R}^{p \times M}$: sparse loading matrix; each column carries exactly $p_1$
   nonzero entries of value $0.9$, drawn without replacement from a shared candidate pool of size
   `overlap_pool_size`.
- $\sigma_m$: amplitude of factor $m$, decreasing over $\{5, 3, 1\}$ (falling back to $1$ for $M > 3$), so
   the components are ordered by strength.
- $\mathcal{S}_1$: PC1's true active loading support; $\widehat{\mathcal{S}}_1$ its estimated counterpart.
- $\boldsymbol{E}$: i.i.d. Gaussian noise, scaled to hit the target SNR.

### Signal-to-noise ratio in decibels

```math
\mathrm{SNR}_{\mathrm{dB}} = 10 \log_{10}\!\left(\frac{\mathrm{Var}(\boldsymbol{Z}\boldsymbol{V}^\top)}{\mathrm{Var}(\boldsymbol{E})}\right).
```

Both the tables and the figures use this decibel axis directly; it is never converted back to a linear ratio.

### False discovery rate, true positive rate, explained variance

```math
\mathrm{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{S}}_1 \setminus \mathcal{S}_1|}{\max\{1, |\widehat{\mathcal{S}}_1|\}}\right],
\qquad
\mathrm{TPR} = \mathbb{E}\left[\frac{|\mathcal{S}_1 \cap \widehat{\mathcal{S}}_1|}{\max\{1, |\mathcal{S}_1|\}}\right].
```

Both are evaluated on the **first PC only**, and FDR control is assumed only there: the factor supports are
 drawn from a shared pool and overlap, so the plug-in PCs that supervise the later components' selectors
 already blend several factors' supports — beyond PC1 there is no unambiguous per-component ground truth to
 control against. (Demo 01's union-FDR section covers what *can* be measured beyond PC1.)

**PEV** — the proportion of explained variance, computed cumulatively across all $M$ components — is the
 quantity the sparsity is bought against, and the place where the cost of FDR control becomes visible.

Classically the explained variance is the plain EV of the estimated components,
 $\mathrm{EV} = \mathrm{tr}(\widehat{\boldsymbol{Z}}^\top \widehat{\boldsymbol{Z}})$. Sparse loading vectors
 are not orthogonal, however, and correlated components double-count shared variance under that trace —
 Zou, Hastie & Tibshirani [[3]](#references) therefore correct it with the **adjusted** EV: a QR
 decomposition of the estimated score matrix, whose diagonal carries each component's *incremental*
 contribution,

```math
\mathrm{EV}_{\mathrm{adj}} = \sum_{m=1}^{M} r_{mm}^2,
\qquad \widehat{\boldsymbol{Z}} = \boldsymbol{Q}\boldsymbol{R}.
```

*What one divides it by* is then what decides what the number means, and this suite reports **both** choices.

**The straightforward reading**, $\mathrm{PEV}_{\mathrm{total}} = \mathrm{EV}_{\mathrm{adj}} /
 \mathrm{Var}(\boldsymbol{X})$, divides by the variance of the whole data matrix. It answers "how much of what
 is in front of us did the $M$ components reproduce", and it is **bounded by 1** by construction — no set of
 $M$ components can explain more variance than $\boldsymbol{X}$ contains.

**The reference paper's reading** (Definition 1 of [[1]](#references)) divides by the **Signal + Mixed EV of
 the data**: the variance carried by the true active variables, with
 $\mathcal{A}$ the union of the factor supports,

```math
\mathrm{PEV}_{\mathrm{sig}} = \frac{\mathrm{EV}_{\mathrm{adj}}}
{\lVert \boldsymbol{X}_{\mathcal{A}} \rVert_F^2 \,/\, (n-1)}.
```

The denominator is **fixed per dataset and shared by all methods** — that is what lets the metric compare
 methods on a common scale, and what the paper's Fig. 3 requires: its curves *rise* with SNR from far below
 $100\%$, which no method-specific (self-normalized) denominator can produce. This convention was
 cross-validated against the published Fig. 3 by the legacy R reference
 (`TRex_Simulations/.../trex_spca/demo_trex_spca_03_fig3.R`, denominator mode `"active"`). An earlier
 revision of the cpp suite normalized each method by its *own* Signal + Mixed EV instead ("Definition 1 as
 printed"), which is structurally $\geq 100\%$ for null-capturing methods and inverts the SNR trend — see
 the [cpp demo-02 README](../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md)
 for the full account.

The consequence is the whole point of the metric: it is **not capped**. A method that spends loadings on null
 variables accumulates variance the denominator does not contain, so its PEV climbs past $100\%$ — and, as
 [[1]](#references) puts it, exceeding $100\%$ "indicates an inferior performance of the respective method".
 Ordinary PCA, using all $p$ loadings, overshoots hardest. Sparse FDR-controlled methods saturate near (but
 below) $100\%$: the gap to $100\%$ is the noise variance sitting on the active columns that $M$ components
 do not capture.

A third reading, $\mathrm{PEV}_{\mathrm{sigmix}}$, reports a method's *own* Signal + Mixed part — its raw EV
 minus the Null EV of its loadings, split by variable — over the same shared denominator. Read for ordinary
 PCA it is the paper's **"Ordinary PCA (Sig + Mix)"** reference curve: it saturates near $100\%$ once all
 components are extracted, while plain ordinary PCA keeps climbing on Null EV.

Keeping all readings is what makes the trade-off legible: the first prices the sparsity, the second exposes
 where a non-sparse method's apparent advantage actually comes from.

---

## T-Rex SPCA specific concepts

- **Plug-in supervision (and the scope of the FDR guarantee)** — each component's selector regresses on
   the *plug-in* ordinary PC $\boldsymbol{z}_m = \boldsymbol{X}\boldsymbol{v}_m$, a one-shot construction
   rather than the iterative refinement most sparse-PCA methods use. The FDR guarantee is therefore
   **conditional on the plug-in PC being a faithful proxy for its factor**: the selector always controls
   the FDR of the regression it is handed, and that transfers to the factor model only where the proxy is
   good. Demo 01 (union-FDR section) and demo 02 (FDR-heatmap section) measure where the condition holds
   and where it fails.
- **`tsn.SPCAMode`** — how the sparse loading vector is assembled once a support has been selected:
  - `ActiveSet`: use only the T-Rex+GVS selected variables' loadings,
  - `Thresholded`: assemble loadings from the full ordinary-PCA solution, restricted to the selected support.
- **`tsn.ENSolverType`** — the elastic-net solver driving each per-component selection: `TENET`
   (Gram-based) or `TENET_AUG` (augmented-LASSO formulation); see
   [../trex_gvs/](../trex_gvs/README.md) for the shared elastic-net background.
- **`lambda_2`** — the ridge penalty (`spca_ctrl.gvs_ctrl.lambda_2`): a negative value selects it by
   $k$-fold cross-validation, $0$ is the degenerate no-ridge case (pure T-LASSO; the wrapper warns), and a
   positive value fixes it.
- **`tFDR`** — the target FDR handed to the per-component T-Rex selector (an argument of
   `sel.select(M, tFDR)`).

---

## Folder contents

```txt
trex_spca/
  ├── README.md
  ├── trex_spca_sim_utils.py           # DGP, MC loops, table + CSV output (Python
  │                                    #   counterpart of trex_spca_sim_utils.hpp)
  ├── trex_spca_plt_utils.py           # shared plotting module (dB sweep figures;
  │                                    #   identical to the cpp suite's module)
  │
  ├── demo_trex_spca_01_mc_sim/
  │   ├── demo_trex_spca_01_mc_sim.py  # the demo source
  │   ├── README.md                    # scenario description and results
  │   ├── generate_plots.sh            # renders this demo's figures from its CSVs
  │   └── simulation_results/
  │       ├── data/                    # .txt summary + .csv table
  │       └── plots/                   # .png + .pdf figures
  │
  └── demo_trex_spca_02_mc_sim_pev/
      ├── demo_trex_spca_02_mc_sim_pev.py
      ├── README.md
      ├── generate_plots.sh            # combines the four sweep CSVs into panels
      └── simulation_results/          # not committed yet — full runs pending
```

The layout mirrors the C++ suite one-to-one.
 [trex_spca_sim_utils.py](trex_spca_sim_utils.py) is the suite-level shared module:
 `generate_sparse_factor_model()` (numpy `default_rng`), `center_columns()`, `evaluate_spca()` (PC1 TPR/FDR
 plus the three PEV readings `PEV`/`PEVsig`/`PEVsigmix`), the parallel MC loops `run_mc_trials_pca()` /
 `run_mc_trials_oracle_pca()` / `run_mc_trials_trex_spca()` (`ProcessPoolExecutor`, default `min(6, cores)`
 workers), and `save_and_print_spca_sweep_results()` (console + `.txt` table with per-sweep-column axis
 precision, plus a tidy pandas CSV `method,metric,<sweep_col>,value`) with the SNR convenience wrapper
 `save_and_print_spca_mc_results()`.

Correctness probes for the individual pipeline stages — TENET vs. TENET_AUG equivalence, scaling
 comparisons, step-by-step C++/R diffs, and a dump of the GVS selector's internal calibration matrices —
 live in the library's test suite under
 `TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`. The C++ and R counterparts of these
 demos run the identical simulation and can be used to cross-check the numbers (on their own RNG streams —
 exact per-trial matches across languages are impossible: numpy PCG64 vs. C++ `std::mt19937`).

---

## Demo descriptions

| # | Name | Purpose | Setting | Main parameters |
|---|------|---------|---------|-----------------|
| **01** | MC Sim | Compare T-Rex SPCA (2 solvers × 2 modes) against OrdPCA/OraclePCA; union-support FDR per PC | Sparse 3-factor model, $n=50$, $p=100$ | $p_1=5$, `overlap_pool_size=30`, tFDR=0.10, SNR sweep in dB (currently one point) |
| **02** | MC Sim PEV | Fig.-3-style cumulative-PEV sweeps (Definition 1 vs. total variance) + FDR/TPR heatmap study | Sparse 3-factor model, $n=50$, $p=100$ | $p_1=10$, sweeps over #PCs / SNR / $p_1$ / tFDR; part 5: (tFDR × SNR) grid |

---

## Running a demo

```bash
cd Python/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim
python demo_trex_spca_01_mc_sim.py                # default min(6, cores) workers
python demo_trex_spca_01_mc_sim.py 8              # 8 worker cores
```

Run inside an environment that carries `trex_selector_neo` (e.g. the `trex-dev` conda environment or the
repo `.venv`). Each demo inserts both its own directory and the parent suite directory onto `sys.path`, so
the shared `trex_spca_sim_utils` module resolves from the nested location and is re-importable by the
spawn-launched Monte Carlo workers — run demos **as files** (not via `python -` / piped stdin), matching the
other Python suites. The optional first CLI argument sets the number of parallel worker cores.

Both demos honour the same environment switches as their cpp counterparts: `TREX_SPCA_NUM_MC=<k>` overrides
the Monte Carlo count for quick pilot runs, and `TREX_SPCA_PARTS` selects the sections/parts to run.

The demos write a human-readable `.txt` summary and tidy `.csv` tables into their local
`simulation_results/data/` folder. Note that only the data is seeded deterministically; each trial's dummies
are drawn from hardware entropy by design, which is what makes the Monte Carlo FDR estimate valid. A re-run
therefore reproduces the reported behaviour to within Monte Carlo noise rather than exactly.

Indices are **0-based** (Python convention); the R counterpart is 1-based.

### Figures

Each demo folder has a `generate_plots.sh` that renders its figures from the saved CSVs into
`simulation_results/plots/` using the shared `trex_spca_plt_utils.py` and the repo-root `.venv`
(matplotlib + pandas + numpy):

```bash
cd demo_trex_spca_01_mc_sim
./generate_plots.sh
```

For demo 01, two figures are rendered: TPR and FDR (against the tFDR target) over the decibel sweep axis,
and — with a `_pev` suffix — the two PEV normalizations side by side, with the 100 % line and OrdPCA's
dashed "Sig + Mix" reference in the Definition-1 panel; a third figure covers the union-FDR CSV when
present. For demo 02, the four sweep CSVs are combined into one Fig.-3-style panel figure per normalization,
plus the part-5 heatmap figure. `--xscale` is accepted for parity with the other suites, but the axis stays
in dB either way.

---

## Notes for new users

- Start with `demo_trex_spca_01_mc_sim/`; move to demo 02 for the explained-variance story.
- Read the PC1-only FDR/TPR rationale above before interpreting results — this is a deliberate
   evaluation-scope choice, not a bug.
- All languages run the per-PC selector with the library default scaling (unit-L2, matching the legacy
   `lm_dummy.R` centering + L2 normalization). Do **not** z-score $\boldsymbol{X}$: the factor amplitude
   signal lives in the column variances, and z-scoring destroys it (measured at $-10$ dB: T-Rex SPCA FDR
   degrades from $\approx 0.14$ to $\approx 0.52$; see the C++
   `validation_trex_spca_06_handrolled_comparison`).
- Cross-check numbers against the language counterparts:
   [`cpp/trex_selector_methods/trex_spca/`](../../../cpp/trex_selector_methods/trex_spca/README.md) and
   `R/trex_selector_methods/trex_spca/`.

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
