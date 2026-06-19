# Simulation results for the `TRex+DA+BT` with equi-correlated hierarchical-block data

## Equi-Correlated Block-Data

The BT route assumes no specific parametric form — it is the **nonparametric case**.
The only assumption is that $\boldsymbol{\Sigma}$ has a hierarchical block structure that can be
well-approximated by a dendrogram.
The DGP that a dendrogram recovers exactly is the **ultrametric covariance model**:
$$
\begin{align*}
[\boldsymbol{\Sigma}]_{jk} &= \rho_{h(j,k)}, \\
h(j,k) &= \text{height of lowest common ancestor of } j \text{ and } k \text{ in the tree}
\end{align*}
$$

In the following we put our interest on a balanced binary tree with two levels and block size $b$, 
this gives the block-Toeplitz structure for

$$
\begin{align*}
\boldsymbol{\Sigma} &= \begin{pmatrix}
\boldsymbol{B}_1 & \rho_2 \mathbf{1} \mathbf{1}^{\top} & \cdots \\
\rho_2 \mathbf{1} \mathbf{1}^{\top} & \boldsymbol{B}_2 & \cdots \\
\vdots & & \ddots \end{pmatrix}\,, \\
\boldsymbol{B}_m &= (1 - \rho_1) \boldsymbol{I}_b + \rho_1 \mathbf{1}_b \mathbf{1}_b^{\top}\,,
\end{align*}
$$

with $\rho_1 > \rho_2$. Valid iff $\rho_1 \in (-\frac{1}{(b-1)},\, 1)$ and $\rho_2 < \rho_1$.
The tree structure is estimated from
$\widehat{\boldsymbol{R}} = \frac{1}{n} \boldsymbol{X}^{\top} \boldsymbol{X}$
via hierarchical clustering of the distance matrix
$$
d_{jk} = 1 - |\hat{\rho}_{jk}|.
$$

$$
\boldsymbol{y} =
\boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon}
$$

$$
\boldsymbol{\varepsilon} \sim \mathcal{N}(0,\, \sigma^2 \boldsymbol{I}_n),
$$
$$
\boldsymbol{\varepsilon} \perp \boldsymbol{X},
$$
with **sparse coefficient vector** $\boldsymbol{\beta} \in \mathbb{R}^p$:
$$
\mathcal{S} = \operatorname{supp}(\beta^*) = \{j : \beta^*_j \neq 0\}, \qquad s = |\mathcal{S}| \ll p
$$
The rows of $X$ are $\text{i.i.d.}$ realisations of a zero-mean, unit-variance Gaussian vector:

$$
\boldsymbol{x}_i = (X_{i,1},\ldots,X_{i,p})^\top \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,\,\boldsymbol{\Sigma}), 
\quad [\boldsymbol{\Sigma}]_{jj} = 1 \;\forall j
$$

---

## SNR Sweep

We examine the scenario for $n=300$, $p=1000$, $s=10$, $\text{blocks}=10$,
$\rho_{\text{within}}=0.5$,
$\rho_{\text{between}}=0.1$,
$200$ Monte Carlos trials per constellation.
The support is generate to represent one active variable per block.
We rely on the `OnePerBlock` policy which generates for the given scenarion the support index as
$\{1, 101, 201, 301, ..., 901 \}$.
For active indices in $\beta_{\text{active}} = 3$, $K=20$ random experiments, and the
target FDR $\alpha := \text{tFDR} = 0.20$.
The $\text{SNR} \in \{0.10, 0.20, 0.50, 1.00, 2.00, 5.00 \}$.
The `TRex+DA+BT` method is used with a `single-linkage` hierarchical clustering algorithm.

|  SNR   | mean_FDP |   mean_TPP |   sd_FDP |     sd_TPP
|--------|----------|------------|----------|-------------
|  0.10  | <span style="background-color: red; color: white;">0.3073</span>   |   0.0050   |   0.4530 |     0.0218
|  0.20  | <span style="background-color: red; color: white;">0.6311</span>   |   0.0255   |   0.4384 |     0.0491
|  0.50  | <span style="background-color: red; color: white;">0.7653</span>   |   0.0880   |   0.2692 |     0.0905
|  1.00  | <span style="background-color: red; color: white;">0.7356</span>   |   0.1535   |   0.2026 |     0.1051
|  2.00  | <span style="background-color: red; color: white;">0.7061</span>   |   0.1925   |   0.1897 |     0.1173
|  5.00  | <span style="background-color: red; color: white;">0.7071</span>   |   0.1965   |   0.1931 |     0.1225

---

## Linkage Sweep

We examine the scenario of $n=300$, $p=1000$, $s=10$, $\text{blocks}=10$,
$\rho_{\text{within}}=0.8$,
$\rho_{\text{between}}=0.1$,
$200$ Monte Carlos trials per constellation.
The support is generate to represent one active variable per block.
We rely on the `OnePerBlock` policy which generates for the given scenarion the support index as
$\{1, 101, 201, 301, ..., 901 \}$.
For active indices in $\beta_{\text{active}} = 3$, $K=20$ random experiments, and the
target FDR $\alpha := \text{tFDR} = 0.20$.
The $\text{SNR} \in \{0.10, 0.20, 0.50, 1.00, 2.00, 5.00 \}$.
The `TRex+DA+BT` method sweeps the hierarchical clustering for three scenarios: `single`,
`complete`, `average`.

### `hc_dist = "single"`

|  SNR   | mean_FDP  |  mean_TPP |   sd_FDP  |    sd_TPP
|--------|-----------|-----------|-----------|------------
|  0.10  | <span style="background-color: red; color: white;">0.2425</span> | 0.0010 | 0.4275 | 0.0100
|  0.20  | <span style="background-color: red; color: white;">0.4887</span> | 0.0075 | 0.4892 | 0.0282
|  0.50  | <span style="background-color: red; color: white;">0.8376</span> | 0.0415 | 0.2783 | 0.0620
|  1.00  | <span style="background-color: red; color: white;">0.8177</span> | 0.0745 | 0.2348 | 0.0814
|  2.00  | <span style="background-color: red; color: white;">0.7663</span> | 0.1135 | 0.2239 | 0.0981
|  5.00  | <span style="background-color: red; color: white;">0.7609</span> | 0.1290 | 0.2010 | 0.0990

### `hc_dist = "complete"`

|  SNR   | ean_FDP |   mean_TPP |   sd_FDP  |    sd_TPP
|--------|---------|------------|-----------|------------
|  0.10  | <span style="background-color: red; color: white;">0.4681</span> | 0.0090 | 0.4871 | 0.0287
|  0.20  | <span style="background-color: red; color: white;">0.8875</span> | 0.0300 | 0.2639 | 0.0521
|  0.50  | <span style="background-color: red; color: white;">0.9259</span> | 0.1390 | 0.0589 | 0.1036
|  1.00  | <span style="background-color: red; color: white;">0.8877</span> | 0.3125 | 0.0548 | 0.1432
|  2.00  | <span style="background-color: red; color: white;">0.8543</span> | 0.5235 | 0.0504 | 0.1601
|  5.00  | <span style="background-color: red; color: white;">0.8413</span> | 0.6935 | 0.0434 | 0.1589

### `hc_dist = "average"`

|  SNR   | mean_FDP |   mean_TPP |   sd_FDP  |    sd_TPP
|--------|----------|------------|-----------|------------
|  0.10  | <span style="background-color: red; color: white;">0.4228</span> | 0.0070 | 0.4834 | 0.0256
|  0.20  | <span style="background-color: red; color: white;">0.8842</span> | 0.0200 | 0.2803 | 0.0425
|  0.50  | <span style="background-color: red; color: white;">0.9215</span> | 0.0970 | 0.0738 | 0.0873
|  1.00  | <span style="background-color: red; color: white;">0.8817</span> | 0.2175 | 0.0718 | 0.1293
|  2.00  | <span style="background-color: red; color: white;">0.8480</span> | 0.3555 | 0.0672 | 0.1565
|  5.00  | <span style="background-color: red; color: white;">0.8358</span> | 0.4635 | 0.0609 | 0.1672

---
