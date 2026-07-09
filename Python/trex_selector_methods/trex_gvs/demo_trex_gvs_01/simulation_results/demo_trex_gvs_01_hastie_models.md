# Simulations `TREX+GVS` - Part 01 (Python)

## Scenarios inspired by testbeds used for the Elastic-Net in Zhu & Hastie (2005)

The following design was proposed in Zhu & Hastie (2005) to evaluate the Elastic Net's ability to
select the entirety of the highly correlated block variables and underline the grouping property
of the elastic network.
We adapt it to evaluate the `TREX+GVS` selector with GVS methods `EN` and `IEN`.

> **Note on indices and numbers.** This Python suite uses **0-based** column indices; the R
> counterpart is 1-based. The numeric result tables below are **reference values** carried over from
> the R / C++ reference implementation (`R/trex_selector_methods/trex_gvs/demo_trex_gvs_01`). Because
> the Python port draws from NumPy's PRNG (`numpy.random.default_rng`) rather than R's RNG, the exact
> per-cell numbers will differ, but the port reproduces the same **qualitative** pattern: `EN` holds
> FDP near `tFDR = 0.1` with TPP rising to 1 as SNR / rho grow, while `IEN` is markedly conservative
> on this all-active design (very low FDP, slowly rising TPP).

We consider the high-dimensional linear regression model:
$$
\boldsymbol{y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon},
$$
$\boldsymbol{y} \in \mathbb{R}^{n}$ is the response vector,
$\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the predictor matrix,
and homoscedastic noise $\boldsymbol{\varepsilon} \in \mathbb{R}^{n}$ with
$\boldsymbol{\varepsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2_{\varepsilon} \boldsymbol{I}_n)$.
The coefficient vector
$\boldsymbol{\beta} \in \mathbb{R}^p$ is sparse true with cardinality
$\text{card}\{\boldsymbol{\beta}\}=s$.

The predictors are structured into $K=3$ active, highly correlated blocks, $\{ G_k \}_{k=1}^{K}$,
each of size $m=50$, driven by independent latent factors.
Using 0-based indexing $j = 0, \ldots, p-1$, the first $p_{\text{active}} = 150$ variables
belong to the active groups, and the remaining $p_{\text{noise}} = p - K \cdot m$ variables are
independent noise.

1. **Latent Factors:** Let $Z_k \sim \mathcal{N}(0, 1)$ for $k \in \{1, 2, 3\}$.
2. **Active Grouped Predictors ($j \le 149$):**
   For each block $G_k$, the predictors are generated as:
   $$
    X_j = Z_k + \xi_j, \qquad \xi_j \sim \mathcal{N}(0, \sigma_x^2),
   $$
   where $\sigma_x^2 = 0.01$. This induces a near-perfect intra-group correlation of
   $\rho = \frac{1}{1 + 0.01} \approx 0.99$.
3. **Null Predictors ($j \ge 150$):**
   $X_j \sim \mathcal{N}(0, 1)$.
4. **True Support ($\beta$):**
   All variables within the three latent groups are genuinely active and contribute equally to the
   response:
   $$
   \beta_j = \begin{cases}
        3 & \text{if } j \le 149, \\
        0 & \text{if } j \ge 150.
    \end{cases}
   $$

The three 0-based active groups are $G_1 = \{0, \ldots, 49\}$, $G_2 = \{50, \ldots, 99\}$,
$G_3 = \{100, \ldots, 149\}$.

## Concepts for the Monte Carlo Demos

### Part 1: Single-Run Validation

* **Objective:** Confirm the integrity of the data generation and demonstrate application of the
  T-Rex+GVS selector.
* **Execution:** Generate one instance of the $200 \times 500$ matrix via `dgp_hastie_snr`.
* **Evaluation:** Run `build_gvs_selector` with `gvs_type` in {`EN`, `IEN`} (and `EN+AUG`).
  Output the single-run TPP and FDP to the console.

---

### Part 2: The Signal-to-Noise Ratio (SNR) Sweep (1D)

* **Objective:** Investigate the robustness under signal-to-noise (SNR) variations.
* **Execution:** Hold the block structure and $\sigma_x^2 = 0.01$ constant.
  We sweep the linear $\text{SNR} \in \{0.1, 0.2, 0.5, 1, 2, 5 \}$ according to
  $$
    \sigma_{\varepsilon} = \sqrt{\frac{\text{Var}(\boldsymbol{X} \boldsymbol{\beta})}{\text{SNR}}}.
  $$
* **Evaluation:** Track `TREX+GVS`'s behaviour as a function of the SNR.

We rely on the empirical signal power $\text{Var}(\boldsymbol{X}\boldsymbol{\beta})$ (with the $n-1$
denominator) to adjust the SNR. For completeness, the theoretical quantities for this scenario are:
for $X_j = Z_k + \xi_j$ we have $\text{Var}[X_j] = 1 + \sigma^2_x$ and, for two members of the same
block, $\text{Cov}[X_i, X_j] = \mathbb{E}[Z_k^2] = 1$. For $\boldsymbol{\beta}_{G_k} = c\,\boldsymbol{1}_m$,
$$
\text{Var}[\boldsymbol{X}_{G_k}\boldsymbol{\beta}_{G_k}]
= c^2\,\boldsymbol{1}^\top_m \boldsymbol{\Sigma}_{G_k} \boldsymbol{1}_m
= c^2\, m\,(m + \sigma_x^2),
$$
so that $\sigma_{\varepsilon} = \sqrt{K c^2 m (m + \sigma_x^2) / \text{SNR}}$.

**Reference — `EN`** ($n=200$, $p=500$, $\text{tFDR}=0.1$, $\sigma_x = 0.1 \Rightarrow \rho \sim 0.99$,
200 Monte Carlo trials per grid point):

| SNR    | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|--------|----------|----------|--------|--------|
| 0.1000 | 0.0705   | 0.4690   | 0.0840 | 0.3641 |
| 0.2000 | 0.1031   | 0.8574   | 0.0787 | 0.2480 |
| 0.5000 | 0.1151   | 0.9984   | 0.0248 | 0.0231 |
| 1.0000 | 0.1176   | 1.0000   | 0.0232 | 0.0000 |
| 2.0000 | 0.1175   | 1.0000   | 0.0247 | 0.0000 |
| 5.0000 | 0.1178   | 1.0000   | 0.0250 | 0.0000 |

**Reference — `IEN`** (same configuration, GVS type `IEN`):

| SNR    | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|--------|----------|----------|--------|--------|
| 0.1000 | 0.0000   | 0.0067   | 0.0000 | 0.0061 |
| 0.2000 | 0.0000   | 0.0126   | 0.0000 | 0.0056 |
| 0.5000 | 0.0000   | 0.0180   | 0.0000 | 0.0035 |
| 1.0000 | 0.0000   | 0.0193   | 0.0000 | 0.0021 |
| 2.0000 | 0.0000   | 0.0198   | 0.0000 | 0.0016 |
| 5.0000 | 0.0000   | 0.0202   | 0.0000 | 0.0012 |

---

### Part 3: The Intra-Group Correlation ($\rho$) Sweep (1D)

* **Objective:** Map the phase transition of the hierarchical clustering due to changes in $\rho$
  for the `TREX+GVS` selector.
* **Execution:** Hold the SNR constant and sweep the within-group noise $\sigma_x^2$ to effectively
  sweep the intra-group correlation $\rho \in [0.1, 0.99]$ (via $\sigma_x = \sqrt{(1-\rho)/\rho}$).
* **Evaluation:** As $\rho$ drops, the blocks become diffuse. We identify the critical $\rho$
  threshold where single-linkage clustering shatters the groups, causing the Elastic Net to
  degenerate into Lasso-like behaviour (a sharp drop in TPP).

**Reference — `EN`** ($n=200$, $p=500$, $\text{tFDR}=0.1$, $\text{SNR}=2.0$, 200 MC per grid point):

| rho  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0673   | 0.4436   | 0.0320 | 0.1042 |
| 0.20 | 0.0841   | 0.8675   | 0.0243 | 0.0597 |
| 0.30 | 0.0913   | 0.9735   | 0.0230 | 0.0261 |
| 0.50 | 0.0942   | 0.9995   | 0.0225 | 0.0029 |
| 0.70 | 0.0946   | 1.0000   | 0.0228 | 0.0000 |
| 0.80 | 0.0944   | 1.0000   | 0.0225 | 0.0000 |
| 0.90 | 0.0939   | 1.0000   | 0.0223 | 0.0000 |
| 0.95 | 0.0942   | 1.0000   | 0.0226 | 0.0000 |
| 0.99 | 0.1174   | 1.0000   | 0.0245 | 0.0000 |

**Reference — `IEN`** (same configuration, GVS type `IEN`):

| rho  | mean_FDP | mean_TPP | sd_FDP | sd_TPP |
|------|----------|----------|--------|--------|
| 0.10 | 0.0435   | 0.3413   | 0.0268 | 0.0798 |
| 0.20 | 0.0302   | 0.7057   | 0.0158 | 0.0730 |
| 0.30 | 0.0287   | 0.8730   | 0.0132 | 0.0503 |
| 0.50 | 0.0231   | 0.9823   | 0.0097 | 0.0175 |
| 0.70 | 0.0229   | 0.9996   | 0.0074 | 0.0015 |
| 0.80 | 0.0258   | 1.0000   | 0.0082 | 0.0000 |
| 0.90 | 0.0292   | 1.0000   | 0.0081 | 0.0000 |
| 0.95 | 0.0307   | 1.0000   | 0.0080 | 0.0000 |
| 0.99 | 0.0000   | 0.0198   | 0.0000 | 0.0016 |

---

### Part 4: 2D Phase Transition Sweep (SNR $\times \rho$)

* **Objective:** Comprehensive picture of the phase transitions of the `TREX+GVS` selector within
  this data scenario.
* **Execution:** Perform a 2D grid search across SNR (rows) and $\rho$ (columns).
* **Evaluation:** Output `mean_TPP` / `mean_FDP` matrices to trace the boundaries separating the
  "safe zone" (FDP $\le 0.10$) from the "failure zone".

#### Elastic Network [EN]

mean_TPP [EN]  (rows: SNR, cols: rho):

|          | rho=0.30 | rho=0.50 | rho=0.70 | rho=0.90 | rho=0.95 | rho=0.99 |
|----------|----------|----------|----------|----------|----------|----------|
| snr=0.20 | 0.164    | 0.478    | 0.697    | 0.821    | 0.830    | 0.858    |
| snr=0.50 | 0.660    | 0.925    | 0.985    | 0.997    | 0.998    | 0.998    |
| snr=1.00 | 0.902    | 0.992    | 0.999    | 1.000    | 1.000    | 1.000    |
| snr=2.00 | 0.974    | 0.999    | 1.000    | 1.000    | 1.000    | 1.000    |
| snr=5.00 | 0.993    | 1.000    | 1.000    | 1.000    | 1.000    | 1.000    |

mean_FDP [EN]  (rows: SNR, cols: rho):

|          | rho=0.30 | rho=0.50 | rho=0.70 | rho=0.90 | rho=0.95 | rho=0.99 |
|----------|----------|----------|----------|----------|----------|----------|
| snr=0.20 | 0.050    | 0.062    | 0.072    | 0.078    | 0.078    | 0.105    |
| snr=0.50 | 0.074    | 0.085    | 0.091    | 0.093    | 0.093    | 0.115    |
| snr=1.00 | 0.084    | 0.093    | 0.095    | 0.095    | 0.094    | 0.117    |
| snr=2.00 | 0.091    | 0.094    | 0.095    | 0.094    | 0.094    | 0.117    |
| snr=5.00 | 0.094    | 0.095    | 0.095    | 0.095    | 0.094    | 0.117    |

#### Informed Elastic Network (IEN)

mean_TPP [IEN] (rows: SNR, cols: rho):

|          | rho=0.30 | rho=0.50 | rho=0.70 | rho=0.90 | rho=0.95 | rho=0.99 |
|----------|----------|----------|----------|----------|----------|----------|
| snr=0.20 | 0.098    | 0.265    | 0.446    | 0.619    | 0.657    | 0.013    |
| snr=0.50 | 0.454    | 0.727    | 0.901    | 0.982    | 0.987    | 0.018    |
| snr=1.00 | 0.721    | 0.923    | 0.989    | 1.000    | 1.000    | 0.019    |
| snr=2.00 | 0.873    | 0.982    | 1.000    | 1.000    | 1.000    | 0.020    |
| snr=5.00 | 0.951    | 0.998    | 1.000    | 1.000    | 1.000    | 0.020    |

mean_FDP [IEN] (rows: SNR, cols: rho):

|          | rho=0.30 | rho=0.50 | rho=0.70 | rho=0.90 | rho=0.95 | rho=0.99 |
|----------|----------|----------|----------|----------|----------|----------|
| snr=0.20 | 0.040    | 0.037    | 0.039    | 0.037    | 0.036    | 0.000    |
| snr=0.50 | 0.046    | 0.047    | 0.044    | 0.044    | 0.046    | 0.000    |
| snr=1.00 | 0.038    | 0.039    | 0.032    | 0.037    | 0.039    | 0.000    |
| snr=2.00 | 0.029    | 0.023    | 0.023    | 0.029    | 0.031    | 0.000    |
| snr=5.00 | 0.017    | 0.012    | 0.017    | 0.020    | 0.021    | 0.000    |

---

*Reference tables adapted from the R/C++ reference implementation. Run the Python demo with the
`RUN_PART_*` flags enabled to reproduce the same qualitative behaviour.*
