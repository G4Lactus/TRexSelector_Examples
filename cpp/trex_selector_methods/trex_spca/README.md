# T-Rex SPCA: Sparse Principal Component Analysis Demonstration Suite

## Overview

This folder introduces the sparse principal component version of the T-Rex selector, the **T-Rex SPCA** according
 to [[1]](#references).
A sparse principal component analysis method built on the T-Rex+GVS elastic-net machinery
 (see [../trex_gvs/](../trex_gvs/README.md)) according to [[2]](#references).

Instead of selecting sparse *variables* in a regression model, T-Rex SPCA selects sparse *loadings* for each
 principal component: the T-Rex+GVS selector is applied per component to control the false discovery rate of
 the estimated loading support. The result is a PCA whose components are supported on a small,
 FDR-controlled set of variables rather than on all $p$ of them.

The demos in this folder are designed to help users understand:

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
 control against.

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
 (`TRex_Simulations/.../trex_spca/demo_trex_spca_03_fig3.R`, denominator mode `"active"`; the R script also
 documents the remaining conventions the paper text leaves open). An earlier revision of this suite
 normalized each method by its *own* Signal + Mixed EV instead ("Definition 1 as printed"), which is
 structurally $\geq 100\%$ for null-capturing methods and inverts the SNR trend — see the demo-02 README
 for the full account.

The goal, in the words of [[1]](#references), is *"to explain the signal and mixed EV with few PCs and sparse
 loadings to allow for interpretability of the obtained PCs"*. Non-sparse PCA methods, or methods that do not
 provide accurate estimates of $\widehat{\boldsymbol{V}}$, are prone to a high Null EV — the variance they
 capture merely corresponds to null (non-active) variables.

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
- **`SPCAMode`** — how the sparse loading vector is assembled once a support has been selected:
  - `ActiveSet`: use only the T-Rex+GVS selected variables' loadings,
  - `Thresholded`: assemble loadings from the full ordinary-PCA solution, restricted to the selected support.
- **`ENSolverType`** — the elastic-net solver driving each per-component selection: `TENET` (Gram-based) or
   `TENET_AUG` (augmented-LASSO formulation); see [../trex_gvs/](../trex_gvs/README.md) for the shared
   elastic-net background.
- **`lambda_2`** — the ridge penalty: a negative value selects it by $k$-fold cross-validation, $0$ disables
   the ridge term, and a positive value fixes it.
- **`tFDR`** — the target FDR handed to the per-component T-Rex selector.

---

## Folder contents

```txt
trex_spca/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex_spca_sim_utils.hpp          # DGP, MC loop, table + CSV output
  ├── trex_spca_plt_utils.py           # shared plotting module (dB sweep figures)
  │
  ├── demo_trex_spca_01_mc_sim/
  │   ├── demo_trex_spca_01_mc_sim.cpp # the demo source
  │   ├── README.md                    # scenario description and results
  │   ├── generate_plots.sh            # renders this demo's figures from its CSV
  │   └── simulation_results/
  │       ├── data/                    # .txt summary + .csv table
  │       └── plots/                   # .png + .pdf figures
  │
  └── demo_trex_spca_02_mc_sim_pev/
      ├── demo_trex_spca_02_mc_sim_pev.cpp
      ├── README.md
      ├── generate_plots.sh            # combines the four sweep CSVs into panels
      └── simulation_results/
          ├── data/                    # one .txt + .csv pair per sweep
          └── plots/                   # .png + .pdf figures
```

Correctness probes for the individual pipeline stages — TENET vs. TENET_AUG equivalence, scaling
 comparisons, step-by-step C++/R diffs, and a dump of the GVS selector's internal calibration matrices —
 live in the library's test suite under
 `TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`, built via its `build_diagnostics`
 target. The R and Python counterparts of this demo run the identical simulation and can be used to
 cross-check the numbers.

---

## Building the demo

From the C++ workspace root:

```bash
cd TRexSelector_Examples/cpp
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/release
```

The executable mirrors the source tree under the `bin/` directory:

```txt
build/release/bin/trex_selector_methods/trex_spca/<demo_folder>/<demo_name>
```

---

## Running a demo

```bash
./build/release/bin/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim
```

The demo writes a human-readable `.txt` summary and a tidy `.csv` table into its local
`simulation_results/data/` folder.

Note that only the data is seeded deterministically; each trial's dummies are drawn from hardware entropy by
 design, which is what makes the Monte Carlo FDR estimate valid. A re-run therefore reproduces the reported
 behaviour to within Monte Carlo noise rather than exactly.

### Figures

The demo folder has a `generate_plots.sh` that renders its figure from the saved CSV into
`simulation_results/plots/` using the shared `trex_spca_plt_utils.py` and the repo-root `.venv`
(matplotlib + pandas + numpy):

```bash
cd demo_trex_spca_01_mc_sim
./generate_plots.sh
```

Two figures are rendered. The first shows TPR and FDR (against the tFDR target) over the decibel sweep axis —
support recovery only. The second — written with a `_pev` suffix — puts the two PEV normalizations side
by side, with the 100 % line drawn in the right panel; it is only produced when the CSV carries the `PEVsig`
metric. When the CSV also carries `PEVsigmix`, OrdPCA's Sig + Mix part is drawn in the Definition-1 panels as
the paper's dashed "Ordinary PCA (Sig + Mix)" reference. `--xscale` is accepted for parity with the other
suites, but the axis stays in dB either way.

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
