# TREX+DA+NN (Nearest Neighbor Dependency Aware T-Rex)

## TREX+DA+NN Model

The dependency-aware relative occurrences are
$$
\Phi^{\text{DA}}_{T,L}(j, \rho_{\text{thr}}) := \Psi_{T,L}(j, \rho_{\text{thr}}) \cdot \Phi_{T,L}(j)
$$
where the penalty factor $\Psi_{T,L}(j, \rho_{\text{thr}}) \in [0.5, 1]$ is defined by
$$
\Psi_{T,L}(j, \rho_{\text{thr}}) := \begin{cases}
\displaystyle \frac{1}{2 - \underset{j' \in \text{Gr}(j, \rho_{\text{thr}} )}{\min} | \Phi_{T,L}(j) - \Phi_{T,L}(j') |} \, , & \quad \text{Gr}(j, \rho_{\text{rho}}) \notin \emptyset \\
\displaystyle 1/2, & \quad \text{Gr}(j, \rho_{\text{thr}}) = \emptyset
\end{cases}
$$
where $\Psi_{T,L}(j, \rho_{\text{thr}})$ penalizes the $j$th ordinary relative occurence according
to its resemblance with other ordinary relative occurrences in its assoicated variable group
$\text{Gr}(j, \rho_{\text{thr}}) \subseteq \{1, \ldots, p \}$.

### The Nearest Neighbors (NN) Penalization

The `TREX+DA+BT` does not allow for arbitrary overlapping groups of highly correlated variables,
which are characteristic for stock return data.

#### Group Design

$$
\text{Gr}(j, \rho_{\text{thr}}) := \{ j' \in \{1, \ldots, p \} \setminus \{j\}:
|\text{corr}(\boldsymbol{x}_{j}, \boldsymbol{x}_{j'})| \geq \rho_{\text{thr}} \}
$$

The challenge is to determine $\rho_{\text{thr}}$ and the other parameters of the dependency-
aware T-Rex selector $( v, T, L )$, such that the `TRex+DA+NN` yields the FDR-controlled solutions
$$
\widehat{\mathcal{A}}_{L}(v, \rho_{\text{thr}}, T) := \{j: \Phi^{NN}_{T,L}(j, \rho_{\text{thr}}) > v \}
$$
where
$$
\Phi^{NN}_{T,L}(j, \rho_{\text{thr}}) := \Psi^{NN}_{T,L}(j, \rho_{\text{thr}}) \cdot \Phi_{T,L}(j)
$$
is the dependency-aware relative occurence of the $j$th variable using the proposed NN group design.

In the following we investigate the `TREX+DA+NN` selector according to the version provided in the
CRAN package `TRexSelector`.

---

## Origin of the term nearest-neighbor

The nearest-neighbor term typically refers to a banded precision matrix (Gaussian Markov field),
not a banded covariance.

The $\text{MA}(\kappa)$ model produces a banded covariance but a dense precision matrix.
Compared to the $\text{AR}(1)$ model the correlation decays faster and truncates sharply at
bandwidth $\kappa$.

On the contrary an $\text{AR}(1)$ process has a tri-diagonal precision matrix but full covariance,
and the correlation dies out slower but still exponential.

---

## Data Setup

### Classical MA$(\kappa)$

A moving-average process of order $\kappa$ expresses an observation as a finite, weighted sum of
white-noise shocks. Let $\{\varepsilon_t\}_{t\in\mathbb{Z}}$ be a white-noise process with
$\mathbb{E}[\varepsilon_t]=0$ and $\operatorname{Var}(\varepsilon_t)=\sigma^2_\varepsilon$.
The causal MA$(\kappa)$ model is
$$
X_t = \mu + \varepsilon_t + \sum_{j=1}^{\kappa} \theta_j \varepsilon_{t-j},
$$
with $\theta_\kappa\neq 0$.
Its covariance function satisfies the hallmark cut-off property $\gamma(k)=0$ for all lags
$k>\kappa$.

### Data-Generation Model

For the simulation study we use a finite-impulse-response filter across the columns of a data
matrix. Draw an $n\times(p+\kappa)$ array of independent innovations
$$
\eta_{i,t}\stackrel{\text{i.i.d.}}{\sim}\mathcal{N}(0,1),
\qquad i=1,\dots,n,\; t=1-\kappa,\dots,p,
$$
and define the $j$-th predictor in the $i$-th sample by the causal convolution
$$
X_{ij}=\sum_{l=0}^{\kappa}\theta_{l}\,\eta_{i,j-l},
\qquad j=1,\dots,p.
\tag{1}
$$

Each row of $\boldsymbol{X}\in\mathbb{R}^{n\times p}$ is therefore an independent realization of
the same $p$-dimensional random vector. In matrix form,
$$
\boldsymbol{X}=\boldsymbol{\eta}\, \boldsymbol{A}^{\mathsf T},
\qquad
\boldsymbol{A}\in\mathbb{R}^{p\times(p+\kappa)},
$$
where $\boldsymbol{A}$ is the banded Toeplitz design matrix
$$
\boldsymbol{A}=
\begin{bmatrix}
\theta_0 &          &        &        &                                    \\
\theta_1 & \theta_0 &        &        &                                    \\
\vdots   & \theta_1 & \ddots &        &                                    \\
\theta_{\kappa} & \vdots   & \ddots & \ddots &                             \\
         & \theta_{\kappa} & \ddots & \ddots & \ddots                      \\
         &          & \ddots & \ddots & \ddots & \ddots                    \\
         &          & & \theta_{\kappa} & \ldots & \theta_{1} & \theta_{0} \\
\end{bmatrix}.
\tag{2}
$$

#### Geometric coefficient parameterisation

To obtain a one-parameter family with controlled overlap, we force the coefficients onto a geometric
path and normalise so that every variable has unit variance for $l=0,\dots,\kappa$
$$
\theta_{l}=c\,\rho^{\,l},
$$
$$
c=\frac{1}{\sqrt{ \left(\sum_{m=0}^{\kappa}\rho^{2m}\right) } }.
\tag{3}
$$

Here $\rho\in\mathbb{R}$ is the **common ratio** of the raw coefficient sequence.
It governs the weighting profile inside the lag window. It is *not* the lag-$1$ autocorrelation,
although for $|\rho|<1$ and $\kappa\gg 1$ the population correlation at band-distance $d$
approaches $\rho^{d}$.

Because the normalisation constant $c$ satisfies
$$
c(\rho)=c(1/\rho)\cdot|\rho|^{-\kappa},
$$
the mapping $\rho\mapsto\{\theta_0,\dots,\theta_\kappa\}$ is essentially symmetric under
$\rho\leftrightarrow 1/\rho$. Consequently the effective correlation-decay regime collapses to
$|\rho|\le 1$; values with $|\rho|>1$ are redundant and are excluded from the simulation.

#### Covariance structure

Let $d=|j-k|$. Since the filter has length $\kappa+1$, variables sharing no innovations are
uncorrelated. For $d\le\kappa$ the covariance is the overlap sum of the coefficient strips:
$$
\Sigma_{jk}
=\sum_{l=0}^{\kappa-d}\theta_{l}\,\theta_{l+d}
=c^{2}\rho^{d}\sum_{l=0}^{\kappa-d}\rho^{2l}.
\tag{4}
$$

Using the geometric-series identity $c^{-2}=\sum_{m=0}^{\kappa}\rho^{2m}$, (4) admits the closed
form
$$
\Sigma_{d}=
\begin{cases}
\displaystyle
\rho^{d}\,
\frac{1-\rho^{2(\kappa-d+1)}}{1-\rho^{2(\kappa+1)}},
& |\rho|\le 1,\;\rho^{2}\neq 1, \\[12pt]
\displaystyle
1-\frac{d}{\kappa+1},
& \rho=+1, \\[12pt]
\displaystyle
(-1)^{d}\Bigl(1-\frac{d}{\kappa+1}\Bigr),
& \rho=-1.
\end{cases}
\tag{5}
$$

The variance is $\Sigma_{0}=1$ by construction. The matrix $\boldsymbol{\Sigma}$ is banded, Toeplitz
(inside the band), and satisfies $\Sigma_{jk}=0$ for $|j-k|>\kappa$.

#### Properties and regimes

| Regime           | Coefficient profile                         | Autocorrelation behaviour |
|------------------|---------------------------------------------|---------------------------|
| $\rho\approx 0$  | $\theta_0\approx 1$, all others negligible. | White-noise-like: correlations vanish immediately. |
| $0<\rho<1$       | Decaying with lag $l$.                      | Roughly geometric decay $\approx\rho^{d}$.         |
| $\rho=1$         | Uniform: $\theta_l=1/\sqrt{\kappa+1}$.      | Triangular (tent-shaped) covariance.               |
| $\rho=-1$        | Alternating sign, uniform magnitude.        | Oscillating triangular covariance.                 |
| $\| \rho \| > 1$ | Redundant: normalisation maps to $\| 1/\rho \| < 1$. | Same as the mirrored $1/\| \rho \|$ case. |

**Weak stationarity** requires no restriction on $\rho$ because the filter is finite.
**Positive definiteness** is immediate from (2):
$\boldsymbol{\Sigma} = \boldsymbol{A} \boldsymbol{A}^{\mathsf \top}$, and
$\theta_0 = c > 0$ guarantees that $\boldsymbol{A}$ has full row rank, so $\boldsymbol{\Sigma}$ is
strictly positive definite for every real $\rho$.

The practically relevant domain for non-oscillatory, locally decaying dependence is
$$
\boxed{|\rho|\le 1},
$$
and most commonly $\rho\in[0,1)$.
The value $\rho=1$ yields a moving-block sum; it is valid but produces the broadest correlation
footprint and the largest condition number for fixed $\kappa$.

---

### Linear Model and SNR Control

We simulate responses from the standard linear model
$$
\boldsymbol{y}=\boldsymbol{X}\boldsymbol{\beta}+\boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon}\sim\mathcal{N}(\mathbf{0},\sigma_{\varepsilon}^{2}\boldsymbol{I}_{n}),
$$
where $\boldsymbol{\varepsilon}$ is the **regression noise**, distinct from the MA innovations $\eta$.
The noise level is set by the signal-to-noise ratio:
$$
\sigma_{\varepsilon}
=\sqrt{\frac{\widehat{\operatorname{Var}}(\boldsymbol{X}\boldsymbol{\beta})}{\text{SNR}}}.
\tag{6}
$$

#### Support Policies

To investigate how banded correlation affects variable selection, active predictors are placed
according to two spatial policies:

1. **CappedSpread.** Fix a maximum gap $\text{max\_gap}$ and set
   $$
   \text{gap}=\min\!\left(\bigl\lfloor p/s \bigr\rfloor,\;\text{max\_gap}\right),
   \qquad
   \mathcal{S}=\bigl\{\,1,\,1+\text{gap},\,1+2\!\cdot\!\text{gap},\,\dots\,\bigr\}.
   $$

    ```R
    gap <- min(p %/% s, max_gap)
    support <- 1L + (0L:(s-1L)) * gap
    ```

2. **Random.** Draw without replacement:
   $$
   \mathcal{S}=\text{sort}\bigl(\text{sample.int}(p,\,\text{size}=s)\bigr).
   $$

    ```R
    sort(sample.int(p, size = s, replace = FALSE))
    ```

---

## SNR Sweep Scenario $\text{MA}(\kappa)$

We examine in simulation `part 2` the performance of the selector for the SNR sweep
$\text{SNR} \in \{0.1, 0.2, \ldots, 1, 2, 5 \}$.
The target FDR $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$,
$\boldsymbol{\beta}_{\text{active}} = 3$, $\rho=0.7$, $200$ Monte Carlo trials,
a Moving-average process order of $\kappa = 3$, and the support is distributed according to
the `CappedSpread(max_gap=20)` policy as
support: {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} over all Monte carlo trials.

| SNR     |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP
|---------|------------|-----------|----------|----------------
| 0.1000  |   0.0000   |  0.0000   |  0.0000  |   0.0000
| 0.2000  |   0.0082   |  0.0010   |  0.0815  |   0.0100
| 0.5000  |   0.1097   |  0.0685   |  0.2427  |   0.1413
| 1.0000  |   <span style="background-color: red; color: white;">0.4528</span>   |  0.3740   |  0.2905  |   0.2265
| 2.0000  |   <span style="background-color: red; color: white;">0.5413</span>   |  0.4600   |  0.2180  |   0.2203
| 5.0000  |   <span style="background-color: red; color: white;">0.5482</span>   |  0.4650   |  0.2087  |   0.2168

For comparison the support generated via `Random` policy per trial:

| SNR     |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP
|---------|------------|-----------|----------|----------------
| 0.1000  |   0.0000   |  0.0005   |  0.0000  |   0.0071
| 0.2000  |   0.0060   |  0.0020   |  0.0661  |   0.0172
| 0.5000  |   0.1124   |  0.0640   |  0.2462  |   0.1319
| 1.0000  |   <span style="background-color: red; color: white;">0.4353</span>   |  0.3140   |  0.3005  |   0.2255
| 2.0000  |   <span style="background-color: red; color: white;">0.5330</span>   |  0.4465   |  0.2376  |   0.2359
| 5.0000  |   <span style="background-color: red; color: white;">0.5510</span>   |  0.4475   |  0.2130  |   0.2246

---

## $\rho$ Sweep Scenario $\text{MA}(\kappa = 3)$

We examine in simulation `part 3` the performance of the selector for the $\rho$ sweep
$\rho \in \{0.1, 0.2, \ldots, 0.9 \}$.
The target FDR $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$,
$\boldsymbol{\beta}_{\text{active}} = 3$, $\text{SNR}=2$, $200$ Monte Carlo trials,
a Moving-average process order of $\kappa = 3$, and the support is distributed according to
the `CappedSpread(max_gap=20)` policy as
support: {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} over all Monte Carlo trials.

| rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP
|------|-----------|-----------|----------|--------------------
| 0.1  |  <span style="background-color: red; color: white;">0.5721</span>   |  0.4295   |  0.1859  |   0.1969
| 0.2  |  <span style="background-color: red; color: white;">0.5806</span>   |  0.4195   |  0.2227  |   0.2119
| 0.3  |  <span style="background-color: red; color: white;">0.6363</span>   |  0.3980   |  0.1735  |   0.1867
| 0.4  |  <span style="background-color: red; color: white;">0.6231</span>   |  0.4100   |  0.1671  |   0.1846
| 0.5  |  <span style="background-color: red; color: white;">0.5809</span>   |  0.4465   |  0.1970  |   0.2086
| 0.6  |  <span style="background-color: red; color: white;">0.5983</span>   |  0.4125   |  0.1934  |   0.1928
| 0.7  |  <span style="background-color: red; color: white;">0.5413</span>   |  0.4600   |  0.2180  |   0.2203
| 0.8  |  0.1648   |  0.8265   |  0.1466  |   0.1499
| 0.9  |  0.1245   |  0.8130   |  0.1177  |   0.1481

For comparison the support generated via `Random` policy per trial:

| rho |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP
|-----|------------|-----------|----------|--------------------
| 0.1 |   <span style="background-color: red; color: white;">0.5780</span>   |  0.4140   |  0.2050  |   0.2143
| 0.2 |   <span style="background-color: red; color: white;">0.5499</span>   |  0.4395   |  0.2112  |   0.2057
| 0.3 |   <span style="background-color: red; color: white;">0.5805</span>   |  0.4405   |  0.1830  |   0.1900
| 0.4 |   <span style="background-color: red; color: white;">0.6020</span>   |  0.4105   |  0.1727  |   0.1789
| 0.5 |   <span style="background-color: red; color: white;">0.5854</span>   |  0.4260   |  0.1941  |   0.2082
| 0.6 |   <span style="background-color: red; color: white;">0.5700</span>   |  0.4195   |  0.2110  |   0.1987
| 0.7 |   <span style="background-color: red; color: white;">0.5330</span>   |  0.4465   |  0.2376  |   0.2359
| 0.8 |   0.1540   |  0.7985   |  0.1440  |   0.1587
| 0.9 |   0.1323   |  0.7690   |  0.1111  |   0.1426

---

## Kappa sweep for $\text{MA}(\kappa)$

At next we examine in simulation `part 4` the performance of the selector for the
$\text{MA}(\kappa)$ scenario.
The target FDR $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$, $\rho=0.7$,
$\boldsymbol{\beta}_{\text{active}} = 3$, $\text{SNR}=2$, $200$ Monte Carlo trials.
The sweep goes over $\kappa \in \{1, 2, 3, 5, 7, 10, 15\}$.

The support is distributed according to the `CappedSpread(max_gap=20)` policy as
support: {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} over all Monte Carlo trials.

| kappa | mean_FDP |  mean_TPP | sd_FDP | sd_TPP 
|-------|----------|-----------|--------|--------
| 1     | <span style="background-color: red; color: white;">0.5938</span> | 0.4285 | 0.1770 | 0.1789
| 2     | <span style="background-color: red; color: white;">0.6100</span> | 0.3945 | 0.1839 | 0.1852
| 3     | <span style="background-color: red; color: white;">0.5213</span> | 0.4765 | 0.2166 | 0.2242
| 5     | <span style="background-color: red; color: white;">0.2810</span> | 0.7115 | 0.2015 | 0.2053
| 7     | <span style="background-color: red; color: white;">0.2165</span> | 0.7815 | 0.1786 | 0.1813
| 10    | <span style="background-color: red; color: white;">0.2102</span> | 0.7850 | 0.1671 | 0.1704
| 15    | <span style="background-color: red; color: white;">0.2313</span> | 0.7925 | 0.1630 | 0.1695

For comparison the support is generated via `Random` policy per trial:

| kappa | mean_FDP |  mean_TPP |  sd_FDP | sd_TPP
|-------|----------|-----------|---------|--------
| 1     | <span style="background-color: red; color: white;">0.5621</span> | 0.4445 | 0.1936 | 0.1864
| 2     | <span style="background-color: red; color: white;">0.5851</span> | 0.3990 | 0.1954 | 0.1926
| 3     | <span style="background-color: red; color: white;">0.5324</span> | 0.4440 | 0.2192 | 0.2216
| 5     | <span style="background-color: red; color: white;">0.2786</span> | 0.7025 | 0.2011 | 0.2070
| 7     | <span style="background-color: red; color: white;">0.2341</span> | 0.7515 | 0.1783 | 0.1936
| 10    | <span style="background-color: red; color: white;">0.2252</span> | 0.7465 | 0.1807 | 0.1883
| 15    | <span style="background-color: red; color: white;">0.2152</span> | 0.7695 | 0.1551 | 0.1819

---

## 2D Sweep: Kappa x Rho

At next we examine in simulation `part 5` the performance of the selector for the 2D sweep of
$\kappa$ vs $\rho$.
The target FDR $\alpha := \text{tFDR} = 0.20$, $n=300$, $p=1000$, $s=10$, $\text{SNR} = 2.0$, 200
MCs.
The parameters:

* $\rho$ sweeps in $\{0.1, 0.3, 0.5, 0.7, 0.8, 0.9\}$
* $\kappa$ sweeps in $\{1, 2, 3, 5, 7, 10, 15\}$

The support is distributed according to the `CappedSpread(max_gap=20)` policy as
support: {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} over all Monte Carlo trials.

 mean_TPP (rows: kappa, cols: rho):

|          |  rho=0.1 | rho=0.3 | rho=0.5 | rho=0.7 | rho=0.8 | rho=0.9 |
|----------|----------|---------|---------|---------|---------|---------|
| kappa=1  |  0.409   | 0.449   | 0.413   | 0.428   | 0.480   | 0.440   |
| kappa=2  |  0.427   | 0.449   | 0.434   | 0.395   | 0.413   | 0.453   |
| kappa=3  |  0.428   | 0.431   | 0.449   | 0.476   | 0.812   | 0.811   |
| kappa=5  |  0.433   | 0.456   | 0.439   | 0.712   | 0.791   | 0.742   |
| kappa=7  |  0.438   | 0.433   | 0.417   | 0.781   | 0.736   | 0.695   |
| kappa=10 |  0.427   | 0.446   | 0.444   | 0.785   | 0.741   | 0.665   |
| kappa=15 |  0.406   | 0.430   | 0.424   | 0.792   | 0.750   | 0.639   |

 mean_FDP (rows: kappa, cols: rho):

|          |  rho=0.1 | rho=0.3 | rho=0.5 | rho=0.7 | rho=0.8 | rho=0.9 |
|----------|----------|---------|---------|---------|---------|---------|
| kappa=1  |  0.592   | 0.582   | 0.613   | 0.594   | 0.544   | 0.584   |
| kappa=2  |  0.572   | 0.585   | 0.588   | 0.610   | 0.574   | 0.532   |
| kappa=3  |  0.568   | 0.596   | 0.573   | 0.521   | <span style="background-color: green; color: white;">0.167</span>   | <span style="background-color: green; color: white;">0.129</span>   |
| kappa=5  |  0.571   | 0.580   | 0.586   | 0.281   | <span style="background-color: green; color: white;">0.155</span>   | <span style="background-color: green; color: white;">0.172</span>   |
| kappa=7  |  0.574   | 0.603   | 0.604   | 0.216   | <span style="background-color: green; color: white;">0.197</span>   | 0.247   |
| kappa=10 |  0.570   | 0.588   | 0.578   | 0.210   | <span style="background-color: green; color: white;">0.203</span>   | 0.272   |
| kappa=15 |  0.599   | 0.607   | 0.602   | 0.231   | 0.224   | 0.389   |

For comparison the support is generated via `Random` policy per trial:

 mean_TPP (rows: kappa, cols: rho):

|          | rho=0.1 | rho=0.3 | rho=0.5 | rho=0.7 | rho=0.8 | rho=0.9 |
|----------|---------|---------|---------|---------|---------|---------|
| kappa=1  | 0.418   | 0.456   | 0.405   | 0.445   | 0.431   | 0.425   |
| kappa=2  | 0.413   | 0.413   | 0.413   | 0.399   | 0.393   | 0.428   |
| kappa=3  | 0.415   | 0.436   | 0.448   | 0.444   | 0.796   | 0.756   |
| kappa=5  | 0.436   | 0.420   | 0.423   | 0.703   | 0.749   | 0.706   |
| kappa=7  | 0.412   | 0.402   | 0.397   | 0.752   | 0.716   | 0.649   |
| kappa=10 | 0.437   | 0.398   | 0.410   | 0.747   | 0.711   | 0.632   |
| kappa=15 | 0.435   | 0.426   | 0.423   | 0.769   | 0.708   | 0.593   |

 mean_FDP (rows: kappa, cols: rho):

|          | rho=0.1 | rho=0.3 | rho=0.5 | rho=0.7 | rho=0.8 | rho=0.9 |
|----------|---------|---------|---------|---------|---------|---------|
| kappa=1  | 0.591   | 0.574   | 0.610   | 0.562   | 0.573   | 0.579   |
| kappa=2  | 0.587   | 0.607   | 0.584   | 0.585   | 0.578   | 0.535   |
| kappa=3  | 0.594   | 0.589   | 0.559   | 0.532   | <span style="background-color: green; color: white;">0.183</span>   | <span style="background-color: green; color: white;">0.149</span>   |
| kappa=5  | 0.571   | 0.612   | 0.594   | 0.279   | <span style="background-color: green; color: white;">0.170</span>   | <span style="background-color: green; color: white;">0.180</span>   |
| kappa=7  | 0.592   | 0.628   | 0.618   | 0.234   | <span style="background-color: green; color: white;">0.183</span>   | 0.253   |
| kappa=10 | 0.566   | 0.627   | 0.602   | 0.225   | 0.221   | 0.300   |
| kappa=15 | 0.577   | 0.602   | 0.589   | 0.215   | 0.250   | 0.396   |

---

## 2D SNR x Rho Sweep Scenarion with $\text{MA}(\kappa = 3)$



---

## New Data Model: AR(1)

In the following we switch the data generating process to an AR(1) model.
For the AR(1) design, the rows of the predictor matrix $\boldsymbol{X} \in \mathbb{R}^{n \times p}$
are generated independently from a $p$-variate normal distribution with mean zero and Toeplitz
covariance matrix $\boldsymbol{\Sigma}$, where
$$
\boldsymbol{\Sigma}_{jk} = \rho^{|j-k|}.
$$
Hence, the correlation between predictors decays geometrically with their index distance.
Given a sparse coefficient vector $\beta$, the response is generated from the linear model
$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
$$
where $\boldsymbol{\varepsilon} \sim \mathcal{N}_n(0, \sigma_\varepsilon^2 \boldsymbol{I}_n)$.
The noise variance is calibrated according to the signal-to-noise ratio,
$$
\sigma_\varepsilon^2 = \operatorname{Var}(\boldsymbol{X}\boldsymbol{\beta})/\mathrm{SNR}
$$
with $\operatorname{Var}(\boldsymbol{X}\boldsymbol{\beta})$ estimated empirically from the
simulated signal.

### 1. Scenario: AR(1) MC — SNR sweep

We exame the performance of the `TREX+DA+NN` selector for $n=300$, $p=1000$, $s=10$, $rho=0.7$,
a target FDR $\alpha := \text{tFDR} = 0.20$, and $\boldsymbol{\beta}_{\text{active}} = 3$.
The support is generated according to the `CappedSpread(max_gap=20)` policy
as {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} for each Monte Carlo trial
and the SNR is swept
over $\text{SNR} \in \{0.1, 0.2, 0.5, 1.0, 2.0, 5.0\}$. For each constellation 200 Monte Carlo
trials are performed.

| SNR    | mean_FDP | mean_TPP | sd_FDP | sd_TPP 
|--------|----------|----------|--------|--------
| 0.1000 | 0.0167   | 0.0010   | 0.1239 | 0.0141 
| 0.2000 | 0.0146   | 0.0055   | 0.1061 | 0.0320 
| 0.5000 | 0.0828   | 0.1145   | 0.1824 | 0.1676 
| 1.0000 | 0.1478   | 0.5730   | 0.1468 | 0.2228 
| 2.0000 | 0.2084   | 0.7980   | 0.1598 | 0.1653 
| 5.0000 | <span style="background-color: red; color: white;">0.2418</span> | 0.7965   | 0.1504 | 0.1676

For comparison the `Random` support policy with redraw per trial is evaluated:

| SNR    | mean_FDP | mean_TPP | sd_FDP | sd_TPP 
|--------|----------|----------|--------|--------
| 0.1000 | 0.0075   | 0.0005   | 0.0789 | 0.0071 
| 0.2000 | 0.0017   | 0.0080   | 0.0236 | 0.0405 
| 0.5000 | 0.1040   | 0.1515   | 0.2053 | 0.1939 
| 1.0000 | 0.1646   | 0.5575   | 0.1532 | 0.2428 
| 2.0000 | <span style="background-color: red; color: white;">0.2101</span>   | 0.7795   | 0.1553 | 0.1797 
| 5.0000 | <span style="background-color: red; color: white;">0.2467</span>   | 0.7725   | 0.1534 | 0.1739 

---

### 2. Scenario: AR(1) MC — Rho sweep

We examine the performance of the `TREX+DA+NN` selector for $n=300$, $p=1000$, $s=10$,
$\text{SNR} = 2.0$,
a target FDR $\alpha := \text{tFDR} = 0.20$, and $\boldsymbol{\beta}_{\text{active}} = 3$.
The support is generated according to the `CappedSpread(max_gap=20)` policy
as {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} for each Monte Carlo trial
and the parameter $\rho$ is swept over $\{0.1, 0.2, ..., 0.9\}$.
For each constellation 200 Monte Carlo trials are performed.

Results — CappedSpread(max_gap=20):

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP 
|-----|----------|----------|--------|--------
| 0.1 | <span style="background-color: red; color: white;">0.5442</span>  | 0.4505  | 0.1925 | 0.1944
| 0.2 | <span style="background-color: red; color: white;">0.5810</span>  | 0.4125  | 0.2317 | 0.2166
| 0.3 | <span style="background-color: red; color: white;">0.6191</span>  | 0.4130  | 0.1761 | 0.1879
| 0.4 | <span style="background-color: red; color: white;">0.6195</span>  | 0.4080  | 0.1874 | 0.1914
| 0.5 | <span style="background-color: red; color: white;">0.5923</span>  | 0.4280  | 0.1934 | 0.2060
| 0.6 | <span style="background-color: red; color: white;">0.5733</span>  | 0.4325  | 0.1921 | 0.1931
| 0.7 | 0.2084   | 0.7980   | 0.1598 | 0.1653 |
| 0.8 | 0.2270   | 0.7365   | 0.1396 | 0.1364 |  
| 0.9 | <span style="background-color: red; color: white;">0.5176</span>   | 0.6030   | 0.1327 | 0.1683

Results — Random support (redrawn per trial):

| rho | mean_FDP | mean_TPP | sd_FDP | sd_TPP 
|-----|----------|----------|--------|--------
| 0.1 | <span style="background-color: red; color: white;">0.5687</span>   | 0.4240   | 0.1996 | 0.2169
| 0.2 | <span style="background-color: red; color: white;">0.5464</span>   | 0.4450   | 0.2228 | 0.2135
| 0.3 | <span style="background-color: red; color: white;">0.6034</span>   | 0.4265   | 0.1813 | 0.1901
| 0.4 | <span style="background-color: red; color: white;">0.6078</span>   | 0.4165   | 0.2008 | 0.2076
| 0.5 | <span style="background-color: red; color: white;">0.6083</span>   | 0.4040   | 0.1879 | 0.1977
| 0.6 | <span style="background-color: red; color: white;">0.5871</span>   | 0.4170   | 0.1926 | 0.1995
| 0.7 | 0.2101   | 0.7795   | 0.1553 | 0.1797 
| 0.8 | 0.2251   | 0.7240   | 0.1376 | 0.1654 
| 0.9 | <span style="background-color: red; color: white;">0.4162</span>   | 0.6050   | 0.1570 | 0.1759

---

### 3. Scenario: AR(1) MC — 2D SNR x Rho sweep

We examine the performance of the `TREX+DA+NN` selector for $n=300$, $p=100$, $s=10$,
$\text{SNR} = 2.0$,
a target FDR $\alpha := \text{tFDR} = 0.20$, and $\boldsymbol{\beta}_{\text{active}} = 3$.
The support is generated according to the `CappedSpread` policy
as {1, 21, 41, 61, 81, 101, 121, 141, 161, 181} for each Monte Carlo trial and the parameter
$\rho$ is swept over $\{0.1, 0.2, \ldots, 0.9\}$ and the $\text{SNR} \in \{0.1, 0.2, 0.5, 1, 2, 5\}$.
For each constellation 200 Monte Carlo trials are evaluated.

Results — CappedSpread(max_gap=20):

 mean_TPP (rows: SNR, cols: rho):

| SNR       | rho=0.1 | rho=0.2 | rho=0.3 | rho=0.4 | rho=0.5 | rho=0.6 | rho=0.7 | rho=0.8 | rho=0.9  
|-----------|---------|---------|---------|---------|---------|---------|---------|---------|----------- 
| SNR=0.10  | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.001   | 0.004   | 0.007    
| SNR=0.20  | 0.000   | 0.001   | 0.000   | 0.000   | 0.000   | 0.000   | 0.005   | 0.032   | 0.047    
| SNR=0.50  | 0.044   | 0.039   | 0.051   | 0.051   | 0.051   | 0.047   | 0.114   | 0.282   | 0.275    
| SNR=1.00  | 0.311   | 0.328   | 0.379   | 0.354   | 0.327   | 0.340   | 0.573   | 0.571   | 0.497    
| SNR=2.00  | 0.451   | 0.412   | 0.413   | 0.408   | 0.428   | 0.432   | 0.798   | 0.737   | 0.603    
| SNR=5.00  | 0.458   | 0.426   | 0.420   | 0.406   | 0.433   | 0.431   | 0.796   | 0.747   | 0.599    

 mean_FDP (rows: SNR, cols: rho):

|          |  rho=0.1 | rho=0.2 | rho=0.3 | rho=0.4 | rho=0.5 | rho=0.6 | rho=0.7 | rho=0.8 | rho=0.9  
|----------|----------|---------|---------|---------|---------|---------|---------|---------|----------- 
| SNR=0.10 |  0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.017   | 0.052   | 0.127    
| SNR=0.20 |  0.000   | 0.002   | 0.000   | 0.005   | 0.005   | 0.000   | 0.015   | 0.100   | 0.350    
| SNR=0.50 |  0.138   | 0.103   | 0.129   | 0.126   | 0.113   | 0.099   | 0.083   | 0.314   | 0.625    
| SNR=1.00 |  <span style="background-color: red; color: white;">0.521</span>   | <span style="background-color: red; color: white;">0.497</span>   | <span style="background-color: red; color: white;">0.522</span>   | <span style="background-color: red; color: white;">0.530</span>   | <span style="background-color: red; color: white;">0.529</span>   | <span style="background-color: red; color: white;">0.489</span>   | 0.148   | <span style="background-color: red; color: white;">0.252</span>   | <span style="background-color: red; color: white;">0.569</span>
| SNR=2.00 |  <span style="background-color: red; color: white;">0.544</span>   | <span style="background-color: red; color: white;">0.581</span>   | <span style="background-color: red; color: white;">0.619</span>   | <span style="background-color: red; color: white;">0.620</span>   | <span style="background-color: red; color: white;">0.592</span>   | <span style="background-color: red; color: white;">0.573</span>   | 0.208   | <span style="background-color: red; color: white;">0.227</span>   | <span style="background-color: red; color: white;">0.518</span>
| SNR=5.00 |  <span style="background-color: red; color: white;">0.546</span>   | <span style="background-color: red; color: white;">0.573</span>   | <span style="background-color: red; color: white;">0.617</span>   | <span style="background-color: red; color: white;">0.628</span>   | <span style="background-color: red; color: white;">0.596</span>   | <span style="background-color: red; color: white;">0.586</span>   | <span style="background-color: red; color: white;">0.242</span>   | <span style="background-color: red; color: white;">0.235</span>   | <span style="background-color: red; color: white;">0.496</span>

For comparison the support is generated according to the `Random` policy and redrawn each trial

mean_TPP (rows: SNR, cols: rho):

|          | rho=0.1 | rho=0.2 | rho=0.3 | rho=0.4 | rho=0.5 | rho=0.6 | rho=0.7 | rho=0.8 | rho=0.9 |
| ---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| SNR=0.10 | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.003   | 0.003   |
| SNR=0.20 | 0.000   | 0.000   | 0.001   | 0.000   | 0.000   | 0.000   | 0.008   | 0.038   | 0.041   |
| SNR=0.50 | 0.029   | 0.051   | 0.067   | 0.066   | 0.054   | 0.062   | 0.151   | 0.280   | 0.233   |
| SNR=1.00 | 0.282   | 0.352   | 0.352   | 0.334   | 0.312   | 0.332   | 0.557   | 0.550   | 0.461   |
| SNR=2.00 | 0.424   | 0.445   | 0.426   | 0.416   | 0.404   | 0.417   | 0.779   | 0.724   | 0.605   |
| SNR=5.00 | 0.435   | 0.446   | 0.424   | 0.420   | 0.414   | 0.425   | 0.772   | 0.729   | 0.612   |

mean_FDP (rows: SNR, cols: rho):

|         | rho=0.1 | rho=0.2 | rho=0.3 | rho=0.4 | rho=0.5 | rho=0.6 | rho=0.7 | rho=0.8 | rho=0.9|
|---------|---------|---------|---------|---------|---------|---------|---------|---------|--------|
|SNR=0.10 | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.000   | 0.008   | 0.047   | 0.121. |
|SNR=0.20 | 0.000   | 0.000   | 0.002   | 0.000   | 0.000   | 0.004   | 0.002   | 0.156   | 0.257. |
|SNR=0.50 | 0.130   | 0.117   | 0.112   | 0.116   | 0.140   | 0.136   | 0.104   | 0.286   | 0.506. |
|SNR=1.00 | <span style="background-color:red; color:white;">0.564</span>   | <span style="background-color:red; color:white;">0.470</span>   | <span style="background-color:red; color:white;">0.521</span>   | <span style="background-color:red; color:white;">0.522</span>   | <span style="background-color:red; color:white;">0.524</span>   | <span style="background-color:red; color:white;">0.475</span>   | 0.165   | <span style="background-color:red; color:white;">0.262</span>    | <span style="background-color:red; color:white;">0.487</span>  |
|SNR=2.00 | <span style="background-color:red; color:white;">0.569</span>   | <span style="background-color:red; color:white;">0.546</span>    | <span style="background-color:red; color:white;">0.603</span>   | <span style="background-color:red; color:white;">0.608</span>    | <span style="background-color:red; color:white;">0.608</span>   | <span style="background-color:red; color:white;">0.587</span>    | 0.210   | <span style="background-color:red; color:white;">0.225</span>    | <span style="background-color:red; color:white;">0.416</span>  |
|SNR=5.00 | <span style="background-color:red; color:white;">0.564</span>   | <span style="background-color:red; color:white;">0.556</span>    | <span style="background-color:red; color:white;">0.610</span>   | <span style="background-color:red; color:white;">0.614</span>    | <span style="background-color:red; color:white;">0.609</span>   | <span style="background-color:red; color:white;">0.584</span>   | <span style="background-color:red; color: white;">0.247</span>   | <span style="background-color:red; color:white;">0.232</span>    | <span style="background-color:red; color:white;">0.392</span> |

---
