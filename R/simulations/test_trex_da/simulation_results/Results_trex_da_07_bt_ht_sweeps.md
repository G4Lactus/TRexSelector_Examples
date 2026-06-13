# Simulation Study: `TREX+DA+BT` with heavy-tailed Block Correlated Data

## Setup

The following simulations examine performance of the `TREX+DA+BT` selector when
$\boldsymbol{X}$ is block-diagonal, consisting of $M$ heavy-tailed correlated blocks driven by the
$p$-dimensional multivariate $t_{\nu}$ distribution with $\nu$ degrees of freedom per block.
The blocks are statistically independent from each other, and we consider the linear model
$$
\boldsymbol{y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{\varepsilon}.
$$
We consider $n$ observations and $p = M \times Q$ variables, within ceteris-paribus simulations
(as far as applicable) to analyze the effects of varying one parameter while fixing all others.
The group correlation structure follows the $Q \times Q$ Toeplitz matrix
$$
\boldsymbol{\Sigma_{m}}(\rho) =
\begin{bmatrix}
1 & \rho & \rho^{2} & \cdots & \rho^{Q-1} \\
\rho & 1 & \rho & \cdots & \rho^{Q-2} \\
\rho^{2} & \rho & 1 & \ldots & \rho^{Q-3} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\rho^{Q-1} & \rho^{Q-2} & \rho^{Q-3} & \cdots & 1
\end{bmatrix},
$$
with
$$
[\boldsymbol{\Sigma}_{m}(\rho)]_{j,k} = \rho^{|j - k|}.
$$

The **total covariance** is then
$$
\boldsymbol{\Sigma} =
\begin{bmatrix}
\boldsymbol{\Sigma}_{1}(\rho) & \boldsymbol{0} & \ldots & \boldsymbol{0} \\
\boldsymbol{0} & \boldsymbol{\Sigma}_{2}(\rho) & \ldots & \boldsymbol{0} \\
\vdots & \vdots & \ddots & \vdots \\
\boldsymbol{0} & \boldsymbol{0} & \ldots & \boldsymbol{\Sigma}_{M}(\rho)
\end{bmatrix}
$$

For the $i$-th row of $\boldsymbol{X}$, we generate column vector
$\boldsymbol{X}_i \in \mathbb{R}^{p}$ by partitioning into $M$ independent blocks
$\boldsymbol{X}_i^{(m)} \in \mathbb{R}^Q$:

$$
\boldsymbol{X}_i =
\begin{pmatrix}
\boldsymbol{X}_i^{(1)} \\ \boldsymbol{X}_i^{(2)} \\ \vdots \\
\boldsymbol{X}_i^{(M)}
\end{pmatrix}.
$$

For the source distribution of $\boldsymbol{X}$ we consider two cases:

* Multivariate Normal distribution: Each block $\boldsymbol{X}^{(m)}_{i}$ uses the AR(1) recursion:
  $$
    z_{j} = \rho \, z_{j-1} + \sqrt{1 - \rho^{2}} \, \varepsilon_{j},
  $$
  for $j = 2, \ldots, Q$ with $z_{1}, \varepsilon_{j} \sim \mathcal{N}(0, 1)$.
  
* Multivariate $t(\nu = 3)$ distribution: Each block transforms via
  $$
   \boldsymbol{X}_{i}^{(m)} = \frac{\boldsymbol{Z}_{i}^{(m)}}{\sqrt{U_{i}^{(m)} / \nu}}
  $$
  where $\boldsymbol{Z}^{(m)}_{i} \sim \mathcal{N}(\boldsymbol{0}, \boldsymbol{\Sigma}_{m}(\rho)))$,
  $U_{i}^{(m)} \sim \chi^2_\nu$, and all $\boldsymbol{Z}_{i}^{(m)}$,
  $U_{i}^{(m')}$ are independent across $m \neq m'$.
  The block covariance scales by $\nu / (\nu - 2)$, yielding the global block-diagonal covariance
  $$
    \text{Cov}(\boldsymbol{X}_i) = \frac{\nu}{\nu - 2} \boldsymbol{\Sigma}.
  $$
  
  **Heavy-Tail Behavior:**
  If a specific block $m$ draws a value of $U_i^{(m)}$ close to zero, all $Q$ variables within that
  specific block will exhibit extreme values simultaneously (tail dependence).
  However, because $U_i^{(1)}, \dots, U_i^{(M)}$ are independent, an
  extreme event in block $m$ has no effect on the magnitude of the variables in the other blocks.

---

## 1. Heavy-tailed X (no white block): SNR Sweep

We consider $n=150$, $p=25$, $M=5$, $Q=5$, $\rho=0.8$, $s=5$ actives,
$\nu=3$ degrees of freedom for the multivariate t-distribution, i.e.,  $t(\nu=3)$,
a target FDR <span style="background-color: blue; color: white">$\alpha := \text{tFDR} = 0.20$</span>, and
$\boldsymbol{\beta}_{\text{active}} = 1.0$, along with $200$ Monte Carlo trials.

### Scenario 1a.: Gaussian Noise

#### 1a.: `hc_dist = "single"`

|  SNR  | mean_FDP         |  mean_TPP | sd_FDP  | sd_TPP  |
|-------|------------------|-----------|---------|---------|
| 0.10  | 0.0425           | 0.0170    | 0.1728  | 0.0627  |
| 0.20  | 0.1283           | 0.0770    | 0.2693  | 0.1483  |
| 0.50  | <span style="background-color: red; color: white;">0.2180</span> | 0.2950 | 0.2944  | 0.2749  |
| 0.60  | <span style="background-color: red; color: white;">0.2152</span> | 0.3690 | 0.2837  | 0.2857  |
| 1.00  | 0.1973           | 0.4950    | 0.2540  | 0.3170  |
| 2.00  | 0.1766           | 0.5620    | 0.2265  | 0.3279  |
| 5.00  | 0.1787           | 0.5540    | 0.2256  | 0.3262  |

#### 1a.: `hc_dist = "complete"`

|  SNR  |  mean_FDP |  mean_TPP |  sd_FDP   |  sd_TPP |
|-------|-----------|-----------|-----------|---------|
|  0.10 |  0.0622   |  0.0220   |  0.2153   |  0.0717 |
|  0.20 |  0.1495   |  0.0990   |  0.2759   |  0.1665 |
|  0.50 |  <span style="background-color: red; color: white;">0.2512</span> | 0.3640 | 0.3001 | 0.2876 |
|  0.60 |  <span style="background-color: red; color: white;">0.2390</span> | 0.4370 | 0.2654 | 0.2836 |
|  1.00 |  <span style="background-color: red; color: white;">0.2346</span> | 0.5890 | 0.2138 | 0.2758 |
|  2.00 |  <span style="background-color: red; color: white;">0.2322</span> | 0.6520 | 0.2129 | 0.2699 |
|  5.00 |  <span style="background-color: red; color: white;">0.2380</span> | 0.6490 | 0.2127 | 0.2686 |

#### 1a.: `hc_dist = "average"`

|  SNR  |  mean_FDP |  mean_TPP |  sd_FDP   |  sd_TPP |
|-------|-----------|-----------|-----------|---------|
|  0.10 |  0.0622   |  0.0220   |  0.2153   |  0.0717 |
|  0.20 |  0.1565   |  0.0930   |  0.2828   |  0.1627 |
|  0.50 |  <span style="background-color: red; color: white;"> 0.2364 </span> |  0.3540 |  0.2892 |  0.2910 |
|  0.60 |  <span style="background-color: red; color: white;"> 0.2387 </span> |  0.4160 |  0.2709 |  0.2929 |
|  1.00 |  <span style="background-color: red; color: white;"> 0.2327 </span> |  0.5650 |  0.2281 |  0.2895 |
|  2.00 |  <span style="background-color: red; color: white;"> 0.2321 </span> |  0.6300 |  0.2201 |  0.2841 |
|  5.00 |  <span style="background-color: red; color: white;"> 0.2316 </span> |  0.6220 |  0.2217 |  0.2876 |

### Scenario 1b.: Heavy $t(\nu=3)$ noise

#### 1b.: `hc_dist = "single"`

|  SNR  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.10 |  0.0820   |  0.0320   |  0.2483  |   0.0991|
|  0.20 |  0.1457   |  0.0970   |  0.2971  |   0.1677|
|  0.50 |  <span style="background-color: red; color: white;">0.2123</span> | 0.3460 | 0.2829 | 0.2965 |
|  0.60 |  <span style="background-color: red; color: white;">0.2129</span> | 0.3800 | 0.2812 | 0.3164 |
|  1.00 |  <span style="background-color: red; color: white;">0.2054</span> | 0.4860 | 0.2727 | 0.3322 |
|  2.00 |  <span style="background-color: red; color: white;">0.2015</span> | 0.5530 | 0.2513 | 0.3264 |
|  5.00 |  0.1912   |  0.5610   |  0.2516  |   0.3299|

#### 1b.: `hc_dist = "complete`

|  SNR  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP |
|-------|-----------|-----------|----------|----------|
|  0.10 |  0.1099   |  0.0430   |  0.2744  |   0.1150 |
|  0.20 |  0.1862   |  0.1230   |  0.3122  |   0.2012 |
|  0.50 |  <span style="background-color: red; color: white;">0.2386</span> | 0.4070 | 0.2622 | 0.2997 |
|  0.60 |  <span style="background-color: red; color: white;">0.2455</span> | 0.4270 | 0.2577 | 0.3045 |
|  1.00 |  <span style="background-color: red; color: white;">0.2567</span> | 0.5570 | 0.2534 | 0.3060 |
|  2.00 |  <span style="background-color: red; color: white;">0.2425</span> | 0.6240 | 0.2259 | 0.2746 |
|  5.00 |  <span style="background-color: red; color: white;">0.2514</span> | 0.6250 | 0.2321 | 0.2835 |

#### 1b.: `hc_dist = "average"`

|  SNR  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP |
|-------|-----------|-----------|----------|----------|
|  0.10 |  0.1086   |  0.0390   |  0.2792  |   0.1093 |
|  0.20 |  0.1885   |  0.1160   |  0.3204  |   0.1874 |
|  0.50 |  <span style="background-color: red; color: white;">0.2417</span>   |  0.3940   |  0.2773  |   0.3047 |
|  0.60 |  <span style="background-color: red; color: white;">0.2495</span>   |  0.4160   |  0.2706  |   0.2990 |
|  1.00 |  <span style="background-color: red; color: white;">0.2495</span>   |  0.5350   |  0.2527  |   0.3106 |
|  2.00 |  <span style="background-color: red; color: white;">0.2396</span>   |  0.5950   |  0.2260  |   0.2867 |
|  5.00 |  <span style="background-color: red; color: white;">0.2342</span>   |  0.6040   |  0.2215  |   0.2884 |

---

## 2. HEAVY-TAILED X (no white block): $\rho$ sweep

We consider $n=150$, $p=25$, $M=5$, $Q=5$, $\rho=$swept, $\text{SNR}=2.0$, $\nu=3.0$,
$\alpha := \text{tFDR} =0.20$, $\boldsymbol{\beta}_{\text{active}}=1.0$, $200$ MC

### Scenario 2a.: Gaussian Noise

#### 2a.: `hc_dist = "single"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0175   |  0.2360   |  0.0838  |   0.2310|
|  0.1  |  0.0121   |  0.2620   |  0.0711  |   0.2307|
|  0.2  |  0.0380   |  0.2930   |  0.1328  |   0.2567|
|  0.3  |  0.0665   |  0.3460   |  0.1841  |   0.2740|
|  0.4  |  0.0894   |  0.4000   |  0.1747  |   0.2926|
|  0.5  |  0.0999   |  0.5280   |  0.1811  |   0.3355|
|  0.6  |  0.1364   |  0.5860   |  0.2065  |   0.3376|
|  0.7  |  0.1366   |  0.6200   |  0.1880  |   0.3276|
|  0.8  |  0.1766   |  0.5620   |  0.2265  |   0.3279|
|  0.9  |  <span style="background-color: red; color: white;">0.2660</span>   |  0.4800   |  0.2780  |   0.3145|

#### 2a.: `hc_dist = "complete"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0822   |  0.6270   |  0.1334  |   0.3341|
|  0.1  |  0.0957   |  0.6340   |  0.1514  |   0.2961|
|  0.2  |  0.1109   |  0.6300   |  0.1643  |   0.3013|
|  0.3  |  0.1201   |  0.6410   |  0.1656  |   0.3108|
|  0.4  |  0.1261   |  0.6670   |  0.1569  |   0.3153|
|  0.5  |  0.1541   |  0.6730   |  0.1966  |   0.3068|
|  0.6  |  0.1566   |  0.6940   |  0.1713  |   0.2849|
|  0.7  |  0.1907   |  0.6760   |  0.1888  |   0.2839|
|  0.8  |  <span style="background-color: red; color: white;">0.2322</span>   |  0.6520   |  0.2129  |   0.2699|
|  0.9  |  <span style="background-color: red; color: white;">0.3193</span>   |  0.5690   |  0.2536  |   0.2934|

#### 2a.: `hc_dist = "average"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0792   |  0.5350   |  0.1485  |   0.3119|
|  0.1  |  0.0678   |  0.5120   |  0.1391  |   0.3168|
|  0.2  |  0.0766   |  0.5360   |  0.1490  |   0.3092|
|  0.3  |  0.1100   |  0.5750   |  0.1580  |   0.3233|
|  0.4  |  0.1294   |  0.5950   |  0.1825  |   0.3192|
|  0.5  |  0.1311   |  0.6560   |  0.1799  |   0.3159|
|  0.6  |  0.1570   |  0.6700   |  0.1899  |   0.3111|
|  0.7  |  0.1832   |  0.6460   |  0.1990  |   0.3072|
|  0.8  |  <span style="background-color: red; color: white;">0.2321</span>   |  0.6300   |  0.2201  |   0.2841|
|  0.9  |  <span style="background-color: red; color: white;">0.3044</span>   |  0.5530   |  0.2548  |   0.2994|

### Scenario 2b.: Heavy-tailed $t(\nu=3)$ Noise

#### 2b.: `hc_dist = "single"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0479   |  0.2250   |  0.1725  |   0.2338|
|  0.1  |  0.0277   |  0.2550   |  0.1178  |   0.2397|
|  0.2  |  0.0414   |  0.2840   |  0.1598  |   0.2586|
|  0.3  |  0.0628   |  0.3330   |  0.1656  |   0.2670|
|  0.4  |  0.0832   |  0.3990   |  0.1810  |   0.3122|
|  0.5  |  0.1112   |  0.5130   |  0.1911  |   0.3441|
|  0.6  |  0.1271   |  0.5720   |  0.1975  |   0.3513|
|  0.7  |  0.1637   |  0.5920   |  0.2353  |   0.3402|
|  0.8  |  <span style="background-color: red; color: white;">0.2015</span>   |  0.5530   |  0.2513  |   0.3264|
|  0.9  |  <span style="background-color: red; color: white;">0.2745</span>   |  0.4690   |  0.2688  |   0.3084|

#### 2b.: `hc_dist = "complete"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0978   |  0.6290   |  0.1683  |   0.3369|
|  0.1  |  0.0800   |  0.6420   |  0.1394  |   0.3116|
|  0.2  |  0.1069   |  0.6290   |  0.1552  |   0.3260|
|  0.3  |  0.1254   |  0.6310   |  0.1923  |   0.3326|
|  0.4  |  0.1278   |  0.6630   |  0.1621  |   0.3212|
|  0.5  |  0.1637   |  0.6600   |  0.1856  |   0.3080|
|  0.6  |  0.1785   |  0.6620   |  0.2017  |   0.2943|
|  0.7  |  <span style="background-color: red; color: white;">0.2325</span>   |  0.6270   |  0.2333  |   0.3019|
|  0.8  |  <span style="background-color: red; color: white;">0.2425</span>   |  0.6240   |  0.2259  |   0.2746|
|  0.9  |  <span style="background-color: red; color: white;">0.3406</span>   |  0.5550   |  0.2493  |   0.2723|

#### 2b.: `hc_dist = "average"`

|  rho  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.0  |  0.0743   |  0.5390   |  0.1405  |   0.3178|
|  0.1  |  0.0705   |  0.5310   |  0.1638  |   0.3149|
|  0.2  |  0.0875   |  0.5310   |  0.1552  |   0.3249|
|  0.3  |  0.1099   |  0.5650   |  0.1926  |   0.3435|
|  0.4  |  0.1228   |  0.5960   |  0.1923  |   0.3325|
|  0.5  |  0.1595   |  0.6400   |  0.2079  |   0.3106|
|  0.6  |  0.1697   |  0.6450   |  0.2084  |   0.3211|
|  0.7  |  <span style="background-color: red; color: white;">0.2187</span>   |  0.6160   |  0.2564  |   0.3266|
|  0.8  |  <span style="background-color: red; color: white;">0.2396</span>   |  0.5950   |  0.2260  |   0.2867|
|  0.9  |  <span style="background-color: red; color: white;">0.3426</span>   |  0.5260   |  0.2595  |   0.2773|

---

## 3. HEAVY-TAILED X (no white block): Q sweep

We consider $n=150$, $p=\text{swept}$, $M=5$, $Q=\text{swept}$, $\rho=0.8$, $\text{SNR}=2.0$,
$\nu=3.0$, $\alpha:=\text{tFDR}=0.20$, $\boldsymbol{\beta}_{\text{active}}=1.0$, $200$ Monte Carlo
trials.

### Scenario 3a.: Gaussian Noise

#### 3a.: `hc_dist = "single"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   0.1766   |  0.5620   |  0.2265  |   0.3279|
|  10 |  50  |   <span style="background-color: red; color: white;">0.2586</span>   |  0.4750   |  0.2606  |   0.3242|
|  15 |  75  |   <span style="background-color: red; color: white;">0.2333</span>   |  0.4940   |  0.2256  |   0.3151|
|  20 |  100 |   <span style="background-color: red; color: white;">0.2887</span>   |  0.4420   |  0.2781  |   0.3230|
|  25 |  125 |   <span style="background-color: red; color: white;">0.2745</span>   |  0.4650   |  0.2728  |   0.3170|
|  30 |  150 |   <span style="background-color: red; color: white;">0.3177</span>   |  0.4550   |  0.2808  |   0.3189|
|  35 |  175 |   <span style="background-color: red; color: white;">0.3673</span>   |  0.4850   |  0.2927  |   0.3128|
|  40 |  200 |   <span style="background-color: red; color: white;">0.3470</span>   |  0.4970   |  0.2805  |   0.3138|
|  45 |  225 |   <span style="background-color: red; color: white;">0.3299</span>   |  0.4840   |  0.2802  |   0.3223|
|  50 |  250 |   <span style="background-color: red; color: white;">0.3941</span>   |  0.4660   |  0.2745  |   0.3048|

#### 3a.: `hc_dist = "complete"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   <span style="background-color: red; color: white;">0.2322</span>   |  0.6520   |  0.2129  |   0.2699|
|  10 |  50  |   <span style="background-color: red; color: white;">0.3051</span>   |  0.5920   |  0.2438  |   0.2980|
|  15 |  75  |   <span style="background-color: red; color: white;">0.3084</span>   |  0.6490   |  0.2227  |   0.2760|
|  20 |  100 |   <span style="background-color: red; color: white;">0.3136</span>   |  0.6050   |  0.2313  |   0.3115|
|  25 |  125 |   <span style="background-color: red; color: white;">0.3464</span>   |  0.6030   |  0.2478  |   0.2774|
|  30 |  150 |   <span style="background-color: red; color: white;">0.3526</span>   |  0.5790   |  0.2454  |   0.3043|
|  35 |  175 |   <span style="background-color: red; color: white;">0.3721</span>   |  0.6270   |  0.2624  |   0.2819|
|  40 |  200 |   <span style="background-color: red; color: white;">0.3587</span>   |  0.6320   |  0.2416  |   0.2803|
|  45 |  225 |   <span style="background-color: red; color: white;">0.4000</span>   |  0.6380   |  0.2537  |   0.2838|
|  50 |  250 |   <span style="background-color: red; color: white;">0.3949</span>   |  0.6350   |  0.2376  |   0.2656|

#### 3a.: `hc_dist = "average"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   <span style="background-color: red; color: white;">0.2321</span>   |  0.6300   |  0.2201  |   0.2841|
|  10 |  50  |   <span style="background-color: red; color: white;">0.2856</span>   |  0.5830   |  0.2409  |   0.2972|
|  15 |  75  |   <span style="background-color: red; color: white;">0.3002</span>   |  0.6250   |  0.2453  |   0.2919|
|  20 |  100 |   <span style="background-color: red; color: white;">0.3120</span>   |  0.5990   |  0.2328  |   0.3044|
|  25 |  125 |   <span style="background-color: red; color: white;">0.3286</span>   |  0.5760   |  0.2405  |   0.2909|
|  30 |  150 |   <span style="background-color: red; color: white;">0.3436</span>   |  0.5630   |  0.2459  |   0.3087|
|  35 |  175 |   <span style="background-color: red; color: white;">0.3696</span>   |  0.6070   |  0.2649  |   0.2894|
|  40 |  200 |   <span style="background-color: red; color: white;">0.3646</span>   |  0.6120   |  0.2408  |   0.2724|
|  45 |  225 |   <span style="background-color: red; color: white;">0.3891</span>   |  0.6120   |  0.2595  |   0.3005|
|  50 |  250 |   <span style="background-color: red; color: white;">0.3890</span>   |  0.6190   |  0.2504  |   0.2783|

### Scenario 3b.: Heavy-tailed $t(\nu=3)$ Noise

#### 3b.: `hc_dist = "single"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   <span style="background-color: red; color: white;">0.2015</span>   |  0.5530   |  0.2513  |   0.3264|
|  10 |  50  |   <span style="background-color: red; color: white;">0.2563</span>   |  0.5120   |  0.2470  |   0.3149|
|  15 |  75  |   <span style="background-color: red; color: white;">0.2542</span>   |  0.4520   |  0.2722  |   0.3386|
|  20 |  100 |   <span style="background-color: red; color: white;">0.2626</span>   |  0.4670   |  0.2585  |   0.3363|
|  25 |  125 |   <span style="background-color: red; color: white;">0.3025</span>   |  0.4410   |  0.3016  |   0.3387|
|  30 |  150 |   <span style="background-color: red; color: white;">0.3161</span>   |  0.4800   |  0.2866  |   0.3252|
|  35 |  175 |   <span style="background-color: red; color: white;">0.3564</span>   |  0.4730   |  0.2900  |   0.3165|
|  40 |  200 |   <span style="background-color: red; color: white;">0.3045</span>   |  0.4830   |  0.2757  |   0.3321|
|  45 |  225 |   <span style="background-color: red; color: white;">0.3414</span>   |  0.4620   |  0.2860  |   0.3134|
|  50 |  250 |   <span style="background-color: red; color: white;">0.3299</span>   |  0.4790   |  0.2890  |   0.3300|

#### 3b.: `hc_dist = "complete"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   <span style="background-color: red; color: white;">0.2425</span>   |  0.6240   |  0.2259  |   0.2746|
|  10 |  50  |   <span style="background-color: red; color: white;">0.3157</span>   |  0.6190   |  0.2111  |   0.2754|
|  15 |  75  |   <span style="background-color: red; color: white;">0.3030</span>   |  0.6300   |  0.2395  |   0.2897|
|  20 |  100 |   <span style="background-color: red; color: white;">0.3055</span>   |  0.5920   |  0.2419  |   0.3118|
|  25 |  125 |   <span style="background-color: red; color: white;">0.3518</span>   |  0.6200   |  0.2505  |   0.2807|
|  30 |  150 |   <span style="background-color: red; color: white;">0.3308</span>   |  0.6180   |  0.2388  |   0.2900|
|  35 |  175 |   <span style="background-color: red; color: white;">0.3538</span>   |  0.6090   |  0.2492  |   0.3036|
|  40 |  200 |   <span style="background-color: red; color: white;">0.3467</span>   |  0.6460   |  0.2427  |   0.2841|
|  45 |  225 |   <span style="background-color: red; color: white;">0.3969</span>   |  0.6310   |  0.2468  |   0.2808|
|  50 |  250 |   <span style="background-color: red; color: white;">0.3939</span>   |  0.6090   |  0.2515  |   0.2901|

#### 3b.: `hc_dist = "average"`

|  Q  |  p   |   mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|------------|-----------|----------|---------|
|  5  |  25  |   <span style="background-color: red; color: white;"> 0.2396 </span>   |  0.5950   |  0.2260  |   0.2867|
|  10 |  50  |   <span style="background-color: red; color: white;"> 0.3131 </span>   |  0.5960   |  0.2339  |   0.2940|
|  15 |  75  |   <span style="background-color: red; color: white;"> 0.2854 </span>   |  0.5880   |  0.2476  |   0.3243|
|  20 |  100 |   <span style="background-color: red; color: white;"> 0.3009 </span>   |  0.5810   |  0.2429  |   0.3058|
|  25 |  125 |   <span style="background-color: red; color: white;"> 0.3559 </span>   |  0.5630   |  0.2749  |   0.3021|
|  30 |  150 |   <span style="background-color: red; color: white;"> 0.3157 </span>   |  0.6090   |  0.2469  |   0.3016|
|  35 |  175 |   <span style="background-color: red; color: white;"> 0.3613 </span>   |  0.6000   |  0.2579  |   0.3119|
|  40 |  200 |   <span style="background-color: red; color: white;"> 0.3415 </span>   |  0.6160   |  0.2574  |   0.3023|
|  45 |  225 |   <span style="background-color: red; color: white;"> 0.4028 </span>   |  0.6020   |  0.2479  |   0.2878|
|  50 |  250 |   <span style="background-color: red; color: white;"> 0.3852 </span>   |  0.5880   |  0.2559  |   0.2998|

---

## 4. HEAVY-TAILED X (no white block): $M$ sweep

We consider $n=150$, $p=\text{swept}$, $Q=5$, $\rho=0.8$, $\text{SNR}=2.0$, $\nu=3.0$,
$\alpha := \text{tFDR}=0.20$, $\boldsymbol{\beta}_{\text{active}}=1.0$, $200$ Monte Carlo trials.

### Scenario 4a.: Gaussian Noise

#### 4a.: `hc_dist = "single"`

|  M  | p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|-----|----|-----------|-----------|----------|---------|
|  1  | 5   | 1  |  0.0675   |  0.2350   |  0.2224  |   0.4251|
|  3  | 15  | 3  |  0.1728   |  0.4150   |  0.2630  |   0.3961|
|  5  | 25  | 5  |  0.1766   |  0.5620   |  0.2265  |   0.3279|
|  10 | 50  | 10 |  <span style="background-color: red; color: white;">0.2589</span> | 0.3970 | 0.2189 | 0.2129 |
|  15 | 75  | 15 |  <span style="background-color: red; color: white;">0.3232</span> | 0.2623 | 0.2336 | 0.1546 |
|  20 | 100 | 20 |  <span style="background-color: red; color: white;">0.3918</span> | 0.1633 | 0.2635 | 0.1133 |
|  30 | 150 | 30 |  <span style="background-color: red; color: white;">0.4164</span> | 0.0963 | 0.2695 | 0.0689 |

#### 4a.: `hc_dist = "complete"`

|  M  | p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|-----|----|-----------|-----------|----------|---------|
|  1  | 5   | 1  |  0.0700   |  0.2400   |  0.2187  |   0.4282|
|  3  | 15  | 3  |  <span style="background-color: red; color: white;">0.2099</span> | 0.5233 | 0.2625 | 0.4072 |
|  5  | 25  | 5  |  <span style="background-color: red; color: white;">0.2322</span> | 0.6520 | 0.2129 | 0.2699 |
|  10 | 50  | 10 |  <span style="background-color: red; color: white;">0.3170</span> | 0.4630 | 0.2088 | 0.2101 |
|  15 | 75  | 15 |  <span style="background-color: red; color: white;">0.3531</span> | 0.2923 | 0.2081 | 0.1482 |
|  20 | 100 | 20 |  <span style="background-color: red; color: white;">0.4441</span> | 0.1953 | 0.2384 | 0.1161 |
|  30 | 150 | 30 |  <span style="background-color: red; color: white;">0.4747</span> | 0.1150 | 0.2360 | 0.0728 |

#### 4a.: `hc_dist = "average"`

|  M  |  p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|----|-----------|-----------|----------|---------|
|  1  |  5   | 1  | 0.0700    |  0.2400   |  0.2187  |   0.4282|
|  3  |  15  | 3  | <span style="background-color: red; color: white;">0.2031</span> |  0.4967   |  0.2575  |   0.4051|
|  5  |  25  | 5  | <span style="background-color: red; color: white;">0.2321</span> |  0.6300   |  0.2201  |   0.2841|
|  10 |  50  | 10 | <span style="background-color: red; color: white;">0.3103</span> |  0.4380   |  0.2203  |   0.2118|
|  15 |  75  | 15 | <span style="background-color: red; color: white;">0.3460</span> |  0.2827   |  0.2070  |   0.1474|
|  20 |  100 | 20 | <span style="background-color: red; color: white;">0.4425</span> |  0.1905   |  0.2496  |   0.1153|
|  30 |  150 | 30 | <span style="background-color: red; color: white;">0.4666</span> |  0.1092   |  0.2482  |   0.0721|

### Scenario 4b.: Heavy-tailed $t(\nu=3)$ Noise

#### 4b.: `hc_dist = "single"`

|  M  |  p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|----|-----------|-----------|----------|---------|
|  1  |  5   | 1  |  0.0425   |  0.2550   |  0.1927  |   0.4370|
|  3  |  15  | 3  |  0.1783   |  0.4017   |  0.2856  |   0.3873|
|  5  |  25  | 5  |  <span style="background-color: red; color:white;">0.2015</span> |  0.5530   |  0.2513  |   0.3264|
|  10 |  50  | 10 |  <span style="background-color: red; color:white;">0.2871</span> |  0.4225   |  0.2272  |   0.2060|
|  15 |  75  | 15 |  <span style="background-color: red; color:white;">0.3298</span> |  0.2593   |  0.2343  |   0.1534|
|  20 |  100 | 20 |  <span style="background-color: red; color:white;">0.4101</span> |  0.1630   |  0.2467  |   0.1081|
|  30 |  150 | 30 |  <span style="background-color: red; color:white;">0.4029</span> |  0.1055   |  0.2626  |   0.0757|

#### 4b.: `hc_dist = "complete"`

|  M  |  p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|----|-----------|-----------|----------|---------|
|  1  |  5   | 1  |  0.0425   |  0.2700   |  0.1792  |   0.4451|
|  3  |  15  | 3  |  <span style="background-color: red; color:white;">0.2100</span> |  0.5300   |  0.2746  |   0.3971|
|  5  |  25  | 5  |  <span style="background-color: red; color:white;">0.2425</span> |  0.6240   |  0.2259  |   0.2746|
|  10 |  50  | 10 |  <span style="background-color: red; color:white;">0.3094</span> |  0.4790   |  0.2007  |   0.2041|
|  15 |  75  | 15 |  <span style="background-color: red; color:white;">0.3680</span> |  0.3090   |  0.2252  |   0.1664|
|  20 |  100 | 20 |  <span style="background-color: red; color:white;">0.4414</span> |  0.2030   |  0.2126  |   0.1199|
|  30 |  150 | 30 |  <span style="background-color: red; color:white;">0.4267</span> |  0.1223   |  0.2516  |   0.0797|

#### 4b.: `hc_dist = "average"`

|  M  |  p   | s  |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-----|------|----|-----------|-----------|----------|---------|
|  1  |  5   | 1  |  0.0425   |  0.2650   |  0.1861  |   0.4424|
|  3  |  15  | 3  |  <span style="background-color: red; color: white;">0.2192</span> |  0.4967   |  0.2885  |   0.3911|
|  5  |  25  | 5  |  <span style="background-color: red; color: white;">0.2396</span> |  0.5950   |  0.2260  |   0.2867|
|  10 |  50  | 10 |  <span style="background-color: red; color: white;">0.3123</span> |  0.4670   |  0.2118  |   0.2033|
|  15 |  75  | 15 |  <span style="background-color: red; color: white;">0.3712</span> |  0.2943   |  0.2182  |   0.1596|
|  20 |  100 | 20 |  <span style="background-color: red; color: white;">0.4373</span> |  0.1958   |  0.2289  |   0.1220|
|  30 |  150 | 30 |  <span style="background-color: red; color: white;">0.4258</span> |  0.1170   |  0.2518  |   0.0749|

---

## 5.  HEAVY-TAILED X (no white block): $\alpha:=\text{tFDR}$ sweep

We consider $n=150$, $p=25$, $M=5$, $Q=5$, $\rho=0.8$, $\text{SNR}=2.0$, $\nu=3.0$,
$\boldsymbol{\beta}_{\text{active}}=1.0$, $200$ Monte Carlo trials.

### Scenario 5a.: Gaussian Noise

#### 5a.: `hc_dist = "single"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color: white;">0.1237</span> |  0.4380   |  0.2128  |   0.3327|
|  0.10 |  <span style="background-color: red; color: white;">0.1466</span> |  0.4890   |  0.2209  |   0.3353|
|  0.15 |  <span style="background-color: red; color: white;">0.1543</span> |  0.5070   |  0.2235  |   0.3378|
|  0.20 |  0.1766   |  0.5620   |  0.2265  |   0.3279|
|  0.25 |  0.2042   |  0.5960   |  0.2357  |   0.3021|
|  0.30 |  0.2083   |  0.6140   |  0.2432  |   0.2881|
|  0.35 |  0.2188   |  0.6290   |  0.2528  |   0.2831|
|  0.40 |  0.2278   |  0.6370   |  0.2551  |   0.2786|
|  0.45 |  0.2328   |  0.6360   |  0.2528  |   0.2703|
|  0.50 |  0.2375   |  0.6370   |  0.2553  |   0.2698|

#### 5a.: `hc_dist = "complete"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color: white;">0.1682</span> |  0.5230   |  0.2100  |   0.3206|
|  0.10 |  <span style="background-color: red; color: white;">0.1830</span> |  0.5880   |  0.2060  |   0.3149|
|  0.15 |  <span style="background-color: red; color: white;">0.2078</span> |  0.5980   |  0.2138  |   0.3087|
|  0.20 |  <span style="background-color: red; color: white;">0.2322</span> |  0.6520   |  0.2129  |   0.2699|
|  0.25 |  <span style="background-color: red; color: white;">0.2524</span> |  0.6760   |  0.2124  |   0.2485|
|  0.30 |  0.2654   |  0.6770   |  0.2170  |   0.2469|
|  0.35 |  0.2757   |  0.6810   |  0.2244  |   0.2432|
|  0.40 |  0.2863   |  0.6780   |  0.2308  |   0.2437|
|  0.45 |  0.3156   |  0.6590   |  0.2444  |   0.2487|
|  0.50 |  0.3207   |  0.6610   |  0.2395  |   0.2474|

#### 5a.: `hc_dist = "average"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color: white;">0.1653</span>   |  0.5180   |  0.2166  |   0.3228|
|  0.10 |  <span style="background-color: red; color: white;">0.1823</span>   |  0.5790   |  0.2049  |   0.3134|
|  0.15 |  <span style="background-color: red; color: white;">0.1958</span>   |  0.5890   |  0.2081  |   0.3146|
|  0.20 |  <span style="background-color: red; color: white;">0.2321</span>   |  0.6300   |  0.2201  |   0.2841|
|  0.25 |  <span style="background-color: red; color: white;">0.2561</span>   |  0.6530   |  0.2165  |   0.2571|
|  0.30 |  0.2642   |  0.6570   |  0.2169  |   0.2426|
|  0.35 |  0.2801   |  0.6600   |  0.2276  |   0.2406|
|  0.40 |  0.2900   |  0.6560   |  0.2368  |   0.2465|
|  0.45 |  0.3126   |  0.6400   |  0.2448  |   0.2528|
|  0.50 |  0.3182   |  0.6440   |  0.2407  |   0.2498|

### Scenario 5b.: Heavy-tailed $t(\nu=3)$ Noise

#### 5b.: `hc_dist = "single"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color:white;">0.1615</span>   |  0.4460   |  0.2573  |   0.3389|
|  0.10 |  <span style="background-color: red; color:white;">0.1762</span>   |  0.4880   |  0.2552  |   0.3371|
|  0.15 |  <span style="background-color: red; color:white;">0.1804</span>   |  0.4960   |  0.2438  |   0.3325|
|  0.20 |  <span style="background-color: red; color:white;">0.2015</span>   |  0.5530   |  0.2513  |   0.3264|
|  0.25 |  0.2123   |  0.5860   |  0.2458  |   0.3011|
|  0.30 |  0.2432   |  0.5990   |  0.2603  |   0.2923|
|  0.35 |  0.2569   |  0.6110   |  0.2566  |   0.2844|
|  0.40 |  0.2735   |  0.6120   |  0.2623  |   0.2833|
|  0.45 |  0.2944   |  0.5980   |  0.2738  |   0.2756|
|  0.50 |  0.2946   |  0.6040   |  0.2684  |   0.2734|

#### 5b.: `hc_dist = "complete"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color:white;">0.1827</span> |  0.5260   |  0.2463  |   0.3229|
|  0.10 |  <span style="background-color: red; color:white;">0.1883</span> |  0.5670   |  0.2344  |   0.3225|
|  0.15 |  <span style="background-color: red; color:white;">0.2072</span> |  0.5720   |  0.2317  |   0.3100|
|  0.20 |  <span style="background-color: red; color:white;">0.2425</span> |  0.6240   |  0.2259  |   0.2746|
|  0.25 |  <span style="background-color: red; color:white;">0.2702</span> |  0.6450   |  0.2337  |   0.2648|
|  0.30 |  <span style="background-color: red; color:white;">0.3077</span> |  0.6340   |  0.2463  |   0.2691|
|  0.35 |  0.3040   |  0.6510   |  0.2359  |   0.2591|
|  0.40 |  0.3190   |  0.6460   |  0.2402  |   0.2635|
|  0.45 |  0.3230   |  0.6440   |  0.2307  |   0.2473|
|  0.50 |  0.3289   |  0.6450   |  0.2332  |   0.2476|

#### 5b.: `hc_dist = "average"`

|  tFDR |  mean_FDP |  mean_TPP |  sd_FDP  |   sd_TPP|
|-------|-----------|-----------|----------|---------|
|  0.05 |  <span style="background-color: red; color:white;">0.1735</span> |  0.5180   |  0.2336  |   0.3178|
|  0.10 |  <span style="background-color: red; color:white;">0.1929</span> |  0.5520   |  0.2342  |   0.3191|
|  0.15 |  <span style="background-color: red; color:white;">0.2117</span> |  0.5480   |  0.2335  |   0.3184|
|  0.20 |  <span style="background-color: red; color:white;">0.2396</span> |  0.5950   |  0.2260  |   0.2867|
|  0.25 |  <span style="background-color: red; color:white;">0.2614</span> |  0.6150   |  0.2342  |   0.2771|
|  0.30 |  <span style="background-color: red; color:white;">0.3014</span> |  0.6100   |  0.2550  |   0.2718|
|  0.35 |  0.312    |  0.6210   |  0.2523  |   0.2715|
|  0.40 |  0.3200   |  0.6180   |  0.2515  |   0.2706|
|  0.45 |  0.3200   |  0.6210   |  0.2426  |   0.2625|
|  0.50 |  0.3264   |  0.6250   |  0.2400  |   0.2621|

---
