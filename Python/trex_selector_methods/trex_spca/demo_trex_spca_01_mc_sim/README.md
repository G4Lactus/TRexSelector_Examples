# Demo 01: T-Rex SPCA on a Sparse Factor Model (Python)

## Purpose

The demo compares **T-Rex SPCA** [[1]](#references) against two PCA baselines on a synthetic sparse
 three-factor model, over a signal-to-noise sweep in decibels.
 Four solver/mode combinations are run — the elastic-net solver `TENET` vs. `TENET_AUG`, each with the
 `ActiveSet` and `Thresholded` loading-assembly modes — against **ordinary PCA** (no sparsity) and
 **oracle-thresholded PCA** (told the true support size).
 The question is whether the T-Rex+GVS machinery controls the false discovery rate of the estimated loading
 support, and what that control costs in explained variance. This demo mirrors **Fig. 2** of
 [[1]](#references) and, like it, uses $p_1 = 5$; the companion
 [demo 02](../demo_trex_spca_02_mc_sim_pev/README.md) mirrors their Fig. 3 at $p_1 = 10$.

Python counterpart of
 [`cpp/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/`](../../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/README.md)
 (same DGP, sweeps, measurement definitions and output-file naming); the R counterpart is
 `R/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/`. Support indices here are **0-based**
 (Python convention).

FDR and TPR are evaluated on **PC1's loading support only**, and FDR *control* is assumed for PC1 only.
 The factor supports are drawn from a shared pool and overlap, and ordinary PCA's orthogonality constraint
 mixes the supports of components 2 and 3 across the true factors — the plug-in PC a later component's
 selector regresses on is already a blend of overlapping supports, so beyond PC1 there is no unambiguous
 per-component ground truth and hence no clean notion of a false discovery to control against. This mirrors
 the reference paper, whose Fig. 2 reports FDR and TPR "for the first PC". PEV does not suffer that
 ambiguity and is cumulative over all $M$ components.

One framing to keep in mind throughout: T-Rex SPCA supervises each component's selector with the *plug-in*
 ordinary PC $\boldsymbol{z}_m = \boldsymbol{X}\boldsymbol{v}_m$ — a one-shot construction, not the
 iterative refinement most sparse-PCA methods use. **The FDR guarantee is therefore conditional on the
 plug-in PC being a faithful proxy for its factor.** Where that proxy is good (PC1, factors above the noise
 floor) control is real; where it degrades, the selector still controls FDR *with respect to the regression
 it was handed*, which is not the same claim. The union-FDR section (Section 2) measures exactly where the
 condition holds.

---

## Data Generation Parameters (`generate_sparse_factor_model`)

We consider a sparse $M$-factor model:

$$
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}
$$

- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the observed data matrix.
- $\boldsymbol{Z} \in \mathbb{R}^{n \times M}$ holds the latent factor scores, column $m$ drawn
   $\mathcal{N}(0, \sigma_m^2)$ with $\sigma_m \in \{5, 3, 1\}$ — the factors are deliberately unequal in
   amplitude.
- $\boldsymbol{V} \in \mathbb{R}^{p \times M}$ is the sparse loading matrix: each column carries exactly
   $p_1$ nonzero entries of value $0.9$.
- $\boldsymbol{E}$ is i.i.d. Gaussian noise, scaled to hit the target SNR in dB.
- $n = 50$, $p = 100$, $p_1 = 5$, $M = 3$. RNG: numpy `default_rng` (PCG64).
- The $p_1$ active indices of each factor are drawn without replacement from a **shared pool** of
   `overlap_pool_size` $= 30$ candidates, so factor supports may partially coincide.

All methods see the same **center-only** $\boldsymbol{X}$ — mean subtraction, no column scaling. That keeps
 ordinary PCA, oracle PCA and T-Rex SPCA on a common covariance-PCA footing. The column scales must *not* be
 normalized here: the factor amplitudes $\sigma_m$ live in the column variances, and z-scoring would destroy
 exactly the signal the factors are distinguished by (see the History note below).

---

## Control Parameters

```text
tFDR = 0.10           # Target FDR for the per-component selector
lambda_2 = -1         # < 0 selects the ridge penalty by k-fold CV (0 = none, > 0 = fixed)
scaling = L2          # Per-component selector scaling mode
MC = 200              # Monte Carlo repetitions per grid point
base_seed = 42        # Trial mc draws its data from base_seed_offset + mc * 1000,
                      # with base_seed_offset = base_seed + round(snr_db)
```

The control object is assembled as `tsn.TRexSPCAControlParameter()` with `mode`
 (`tsn.SPCAMode.ActiveSet` / `Thresholded`), `en_solver` (`tsn.ENSolverType.TENET` / `TENET_AUG`) and the
 nested `gvs_ctrl` (`lambda2_method = tsn.LambdaSelectionMethod.CV_1SE_CCD`, `lambda_2 = -1.0`,
 `trex_ctrl.scaling_mode = tsn.ScalingMode.L2`). The selector is the root-level wrapper
 `tsn.TRexSPCASelector(X, spca_ctrl=ctrl, seed=-1)`; `sel.select(M, tFDR)` returns `self`, and the results
 are read from the `sel.V`, `sel.Z` and `sel.active_sets` properties (0-based). The PCA baselines use
 `trex_selector_neo.ml_methods.PCA(center=False, normalize=False).fit(X, M)` → `PCAResult` with `.V` / `.Z`
 — `center=False` because $\boldsymbol{X}$ is already centered (covariance PCA).

Note that only the **data** is seeded deterministically. Each trial's dummies are drawn from hardware entropy
 (selector seed $-1$) by design, which is what makes the Monte Carlo FDR estimate valid — so a re-run
 reproduces the committed numbers only to within Monte Carlo noise ($\approx \pm 0.01$ at 200 trials), not
 exactly.

---

## Methods Compared

| Method | Sparsity mechanism | `SPCAMode` | `ENSolverType` |
| --- | --- | --- | --- |
| **OrdPCA** | none — all $p$ loadings retained | — | — |
| **OraclePCA** | top-$p_1$ by magnitude, true support size known | — | — |
| **TRexSPCA-EN-Act** | T-Rex+GVS selection [[2]](#references) | `ActiveSet` | `TENET` |
| **TRexSPCA-ENaug-Act** | T-Rex+GVS selection [[2]](#references) | `ActiveSet` | `TENET_AUG` |
| **TRexSPCA-EN-Thr** | T-Rex+GVS selection [[2]](#references) | `Thresholded` | `TENET` |
| **TRexSPCA-ENaug-Thr** | T-Rex+GVS selection [[2]](#references) | `Thresholded` | `TENET_AUG` |

The two modes differ in how the loadings are rebuilt once the support is chosen: `ActiveSet` keeps the
 selector's ridge coefficients, while `Thresholded` re-solves the ordinary PCA problem restricted to the
 selected support.

---

## Metrics

FDR and TPR are evaluated on PC1's loading support. With $\mathcal{S}_1$ the true active support of
 factor 1 and $\widehat{\mathcal{S}}_1$ the estimated (nonzero) support of the first sparse loading vector,

```math
\mathrm{FDR} = \mathbb{E}\left[\frac{|\widehat{\mathcal{S}}_1 \setminus \mathcal{S}_1|}
{\max\{1, |\widehat{\mathcal{S}}_1|\}}\right],
\qquad
\mathrm{TPR} = \mathbb{E}\left[\frac{|\mathcal{S}_1 \cap \widehat{\mathcal{S}}_1|}
{\max\{1, |\mathcal{S}_1|\}}\right],
```

with the expectations estimated as Monte Carlo means over the trials.

For the explained variance, the classical quantity is the EV of the estimated components,
 $\mathrm{EV} = \mathrm{tr}(\widehat{\boldsymbol{Z}}^\top \widehat{\boldsymbol{Z}})$. Sparse PCA drops the
 orthogonality constraint, so correlated components double-count shared variance under this trace.
 Zou, Hastie & Tibshirani [[3]](#references) therefore correct it: with
 $\widehat{\boldsymbol{Z}} = \boldsymbol{Q}\boldsymbol{R}$ the QR decomposition of the estimated score
 matrix, the **adjusted EV** keeps only each component's incremental contribution,
 $\mathrm{EV}_{\mathrm{adj}} = \sum_{m=1}^{M} r_{mm}^2$.

The **TRex-SPCA** paper [[1]](#references) adopts this corrected EV and builds its **Definition 1** on top
 of it, asking *where* the variance comes from: variance on the true active variables (Signal + Mixed EV) is
 worth explaining, variance on null variables (Null EV) is not. Accordingly, the demo reports **three**
 readings:

- The share of the **total** variance, bounded by 1:
  $\mathrm{PEV} = \mathrm{EV}_{\mathrm{adj}} / \mathrm{Var}(\boldsymbol{X})$.
- **Definition 1** of [[1]](#references):
  $\mathrm{PEV}_{\mathrm{sig}} = \mathrm{EV}_{\mathrm{adj}} / \bigl(\lVert \boldsymbol{X}_{\mathcal{A}}
   \rVert_F^2 / (n-1)\bigr)$, where the denominator is the **Signal + Mixed EV of the data** — the variance
   carried by the true active variables ($\mathcal{A}$ = union of the factor supports), fixed per dataset
   and **shared by all methods**. Methods that additionally capture null-variable variance push past
   $100\%$, so this ratio is **not** capped at 1 — exceeding $100\%$ is the failure the metric is designed
   to expose. (Denominator convention validated against the published Fig. 3 by the legacy R reference,
   `demo_trex_spca_03_fig3.R`, mode `"active"`; it is *not* the method's own Signal + Mixed EV, which is
   structurally $\geq 100\%$ for null-capturing methods and inverts the SNR trend — see the
   [cpp demo-02 README](../../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md)
   for the full account.)
- $\mathrm{PEV}_{\mathrm{sigmix}}$ — the method's *own* Signal + Mixed part over the same shared
   denominator; read for OrdPCA it is the paper's "Ordinary PCA (Sig + Mix)" reference curve.

The CSV tags them `PEV`, `PEVsig` and `PEVsigmix`.

---

## The Sweep

The cpp reference runs a single **SNR sweep in decibels** over the legacy CRAN grid
 $\mathrm{SNR} \in \{-10, -7, -5, -3, 0, 3, 5, 7, 10\}$ dB, 200 MC trials per point.
 **This Python port currently runs the single point $-10$ dB** (`snr_values = [-10.0]` in `main()`, a
 deliberate fast-validation downscale — one grid point of the full sweep already takes several minutes per
 method in Python); the committed results match this single-point configuration. Extending `snr_values` to
 the full grid makes the demo identical to the cpp reference.

Section 2 reruns the same grid for the **union-support FDR per PC** (`ActiveSet`/`TENET` only) — see the
 [dedicated section below](#union-support-fdr-beyond-pc1) for what it measures and why the union of the
 factor supports is the only well-posed per-PC truth beyond PC1.

---

## Output Files

Written to `simulation_results/data/`:

- `demo_trex_spca_01_mc_sim.txt` / `.csv` — FDR, TPR and the PEV readings per method and SNR level. The CSV
   is tidy long format `method,metric,snr_db,value` with metrics `FDR`, `TPR`, `PEV` (total-variance
   normalization), `PEVsig` (Definition 1) and `PEVsigmix` (the method's Sig + Mix part).
- `demo_trex_spca_01_mc_sim_union_fdr.csv` — Section 2: union-support FDR and mean support size per PC
   (`pc,metric,snr_db,value`, metrics `FDRunion` / `k`).

Figures (PNG + PDF) go to `simulation_results/plots/`, produced by `./generate_plots.sh` via the shared
 suite plotting module [`../trex_spca_plt_utils.py`](../trex_spca_plt_utils.py):

- `demo_trex_spca_01_mc_sim.png` — TPR and FDR over the sweep (support recovery only).
- `demo_trex_spca_01_mc_sim_union_fdr.png` — Section 2: union-support FDR per PC.
- `demo_trex_spca_01_mc_sim_pev.png` — the two PEV normalizations side by side (the explained-variance
   story lives entirely in this figure).

---

## Running the Demo

```bash
cd Python/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim
python demo_trex_spca_01_mc_sim.py                # default min(6, cores) workers
python demo_trex_spca_01_mc_sim.py 8              # 8 worker cores
./generate_plots.sh                               # render the figures from the saved CSVs
```

`TREX_SPCA_NUM_MC=<k>` overrides the Monte Carlo count for quick pilot runs (e.g. `TREX_SPCA_NUM_MC=10`),
and `TREX_SPCA_PARTS` selects the sections (default `12`; `1` = method comparison, `2` = union-FDR sweep) —
the same environment switches as the cpp executable.

Run the demo inside an environment that carries `trex_selector_neo` (e.g. the `trex-dev` conda environment
or the repo `.venv`). The Monte Carlo trials run in parallel via `concurrent.futures.ProcessPoolExecutor`
(spawn), so the demo must be run **as a real `.py` file** (workers re-import the module) — not via
`python -` / piped stdin. The optional first CLI argument sets the number of worker cores.

---

## Union-support FDR beyond PC1

What can be said about the later components?
A per-factor FDR cannot be controlled!
It would brand legitimately selected leaked variables as false discoveries and measure the ambiguity of
 the truth rather than the selector.
The **union of the factor supports** is well-posed for every component: a variable outside
 the union is null for *every* factor, so selecting it is a false discovery no matter how the components
 mix.
Section 2 of the demo (`TREX_SPCA_PARTS=2`) therefore scores each PC's selected support
 (`ActiveSet`/`TENET`, the selector's own active sets) against that union — answering the one question the
 overlap leaves answerable: *does the selector start harvesting pure-noise variables on the later PCs?*

Formally, with $\mathcal{S}_{\cup} = \bigcup_{k=1}^{M} \mathcal{S}_k$ the union of the true factor supports
 and $\widehat{\mathcal{A}}_m$ the selector's active set for component $m$,

```math
\mathrm{FDR}_{\cup}(m) = \mathbb{E}\left[\frac{|\widehat{\mathcal{A}}_m \setminus \mathcal{S}_{\cup}|}
{\max\{1, |\widehat{\mathcal{A}}_m|\}}\right],
```

again as a Monte Carlo mean; the CSV also records the companion support size
 $\mathbb{E}\bigl[|\widehat{\mathcal{A}}_m|\bigr]$. For $m = 1$ this differs from the FDR above only in the
 more lenient truth ($\mathcal{S}_1 \subseteq \mathcal{S}_{\cup}$), which is why PC1's union-FDR sits
 slightly below its per-factor FDR.

The cpp reference (full grid, 200 MC) finds: **PC1 is controlled everywhere**; **PC2 and PC3 are controlled
 exactly where their factors are alive, and break where they drown** — the factors' amplitude offsets
 ($\approx 4$ dB and $\approx 14$ dB below factor 1) predict precisely where each PC's curve breaks. Once a
 factor is buried, its plug-in PC is essentially a noise direction, and the selector — asked to explain a
 noise response — dutifully selects noise variables, with a correspondingly inflated support size. The
 Python smoke checks reproduce this pattern at $-10$ dB (PC1 controlled, PC2/PC3 broken); see the
 [cpp demo-01 README](../../../../cpp/trex_selector_methods/trex_spca/demo_trex_spca_01_mc_sim/README.md)
 for the quantitative account. **Consequence**: per-component control is only as good as the component's
 effective SNR — a caveat inherited from the one-shot plug-in construction, not from the selector.

---

## Committed Results (SNR = -10 dB, num_MC = 200, run 2026-07-08)

| Method | FDR | TPR | PEV |
|---|---|---|---|
| OrdPCA | 0.9500 | 1.0000 | 0.2014 |
| OraclePCA | 0.0080 | 0.9920 | 0.1233 |
| TRexSPCA-EN-Act | 0.1317 | 0.9970 | 0.1202 |
| TRexSPCA-ENaug-Act | 0.1350 | 0.9970 | 0.1203 |
| TRexSPCA-EN-Thr | 0.1359 | 0.9980 | 0.1335 |
| TRexSPCA-ENaug-Thr | 0.1354 | 0.9980 | 0.1331 |

> **Note**: this table predates the 2026-07-21 rework — it was produced before the demo reported the
> `PEVsig` / `PEVsigmix` readings and the union-FDR section, so the committed CSV carries only
> `FDR`/`TPR`/`PEV` and there is no committed union-FDR CSV yet. The reworked code reproduces the same
> FDR/TPR/PEV values to within MC noise (verified by pilot runs); a full re-run to refresh the committed
> files with the new columns is pending.

### Interpretation

- **OrdPCA**'s $\mathrm{FDR}\approx0.95$ and $\mathrm{TPR}=1.0$ are trivial artifacts of retaining all
   $p=100$ loadings (true support is $p_1=5$) — a non-sparse reference line, not a competitor.
- **OraclePCA** achieves near-perfect PC1 recovery (it is told the true support size); its PEV ($0.123$) is
   the reference for how much variance the true sparse signal explains, independent of any selection
   procedure.
- All four **T-Rex SPCA** variants recover PC1's support almost perfectly (TPR $\approx 0.997$–$0.998$)
   with FDR close to the $\mathrm{tFDR}=0.10$ target — even though $-10$ dB sounds like a very weak signal,
   the strong first factor (std $5$) dominates PC1 recovery on the covariance-PCA footing.
- The realized FDR ($\approx 0.13$–$0.14$) sits $\approx 0.03$ above the target and the legacy CRAN R
   reference ($0.100$). This is the long-known residual gap between the C++ core and the CRAN R
   implementation, investigated by the cross-check programs in the TRexSelector library test suite
   (`TRexSelector/cpp/tests/validation/trex_selector_methods/trex_spca/`);
   `validation_trex_spca_06_handrolled_comparison` confirms it is **not** a `TRexSPCASelector` artifact.
- The four T-Rex variants are statistically indistinguishable from each other at this single SNR point
   (the cpp full-grid run confirms `TENET` and `TENET_AUG` are the same estimator on this problem — given
   identical dummy seeds they produce bit-identical selections). The **Thresholded** variants' slightly
   higher PEV ($\approx 0.133$ vs. $\approx 0.120$) is the one systematic difference: re-solving on the
   selected support recovers more variance than the active-set ridge coefficients do.

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline
> z-score-standardized $\boldsymbol{X}$ and set the per-PC selector to `ScalingMode.ZSCORE`, which converts
> the covariance PCA into a correlation PCA and destroys the factor amplitude signal — the measured T-Rex
> FDR inflated to $\approx 0.52$ at this same nominal SNR. The C++ validation program
> `validation_trex_spca_06_handrolled_comparison` isolated the cause, and the pipeline is back to
> center-only.

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
