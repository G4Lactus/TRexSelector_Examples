# Classical T-Rex Selector: Demonstration Suite

## Overview

This folder contains C++ example programs for the classical **Terminating Random Experiments Selector (T-Rex Selector)**.

The demos are designed to help first-time users understand three things:

1. how to generate synthetic regression data,
2. how to run the classical T-Rex Selector in different settings,
3. how to interpret the resulting variable selections and simulation summaries.

If you are new to this folder, start with **Demo 01** for a basic single run and then continue with **Demo 02** for a Monte Carlo study across several signal-to-noise ratios.

---

## What this folder covers

The examples study sparse variable selection in Gaussian linear regression in both:

- **low-dimensional settings**, where $(n > p)$,
- **high-dimensional settings**, where $(p > n)$.

Across the demos, we vary:

- the sparsity pattern of the true coefficient vector,
- the signal strength,
- the support construction,
- and the computational strategy (standard in-memory or memory-mapped workflows).

The main statistical goal is to recover the active variables while controlling the false discovery rate.

---

## Statistical model and notation

We work with the classical Gaussian linear model

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon},
\qquad
\boldsymbol{\varepsilon} \sim \mathcal{N}(0, \sigma^2 \boldsymbol{I}_n).
$$

### Notation

- $n$: Number of observations.
- $p$: Number of variables.
- $\boldsymbol{y} \in \mathbb{R}^n$: Response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$: Design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$: Regression coefficient vector, $\beta_{j}$ denotes the $j$-th coefficient.
- $\boldsymbol{\varepsilon} \in \mathbb{R}^n$: Gaussian noise vector.
- $s = |\mathcal{A}|$: Number of truly active variables.
- $\mathcal{A} = \{j : \beta_j \neq 0\}$: True support set.
- $\widehat{\mathcal{A}}$: Estimated support set returned by the selector.

### Data-generating variations

Depending on the demo, the simulation setup may use:

- **fixed or random support indices**,
- **deterministic nonzero coefficients** such as $\beta_j = 1$,
- **random nonzero coefficients**, for example drawn from $\mathcal{N}(0,1)$.

### Signal-to-noise ratio

The linear signal-to-noise ratio is defined as

$$
\mathrm{SNR} = \frac{\mathrm{Var}(\boldsymbol{X}\boldsymbol{\beta})}{\mathrm{Var}(\boldsymbol{\varepsilon})}.
$$

For a fixed design matrix $\boldsymbol{X}$ and coefficient vector $\boldsymbol{\beta}$, the noise variance $\sigma^2$ is adjusted to match a desired target SNR.

Typical SNR values used in the demos are

$$
\{0.1,\; 0.5,\; 1,\; 2,\; 5\},
$$

which range from weak to strong signal settings.

---

## Statistical targets

The selector aims to identify relevant variables while controlling false discoveries.

### False Discovery Rate (FDR)

$$
\mathrm{FDR}
=
\mathbb{E}
\left[
\frac{|\widehat{\mathcal{A}} \setminus \mathcal{A}|}
{\max\{1, |\widehat{\mathcal{A}}|\}}
\right].
$$

In the simulations, the target level is denoted by **tFDR**, typically

$$
\alpha = 0.1.
$$

### True Positive Rate (TPR)

$$
\mathrm{TPR}
=
\mathbb{E}
\left[
\frac{|\mathcal{A} \cap \widehat{\mathcal{A}}|}
{\max\{1, |\mathcal{A}|\}}
\right].
$$

Informally:

- **FDR** measures how many selected variables are false discoveries on average.
- **TPR** measures how many truly active variables are recovered on average.

A good procedure keeps FDR below the target level while achieving high TPR.

---

## Start here

For first use, the following order is recommended:

1. **Demo 01 — Single Run**  
   Best starting point. Shows the basic workflow in one low-dimensional and one high-dimensional example.

2. **Demo 02 — Monte Carlo, Fixed Support**  
   Shows how performance changes across SNR values and provides aggregated FDR/TPR output.

3. **Demo 05 — Memory-Mapped Single Run**  
   Useful if you are interested in large-data workflows or limited-RAM settings. To reveal its full power, we recommend using an SSD drive.

---

## Demonstration structure

Each demo subfolder typically contains:

- **`demo_*.cpp`**: the source file for the demo.
- **`README.md`**: a description of the scenario and the main interpretation.
- **`simulation_results/`**: generated output files such as text summaries and CSV tables.

This makes each demo self-contained and easier to inspect independently.

---

## Demo descriptions

| # | Name | Purpose | Setting | Main parameters |
|---|------|---------|---------|-----------------|
| **01** | Single Run | Demonstrates basic T-Rex usage | One low-dimensional and one high-dimensional scenario | $n \in \{150, 5000\}$, $p \in \{300, 1000\}$, $s = 6$, SNR = 1.0, tFDR = 0.1 |
| **02** | MC: Fixed Support | Studies FDR/TPR over an SNR range | Moderate dimension, fixed support | $n = 300$, $p = 1000$, $s = 10$, SNR $ \in \{0.1, 0.5, 1, 2, 5\} $, tFDR = 0.1, MC = 200 |
| **03** | MC: Variable Support | Studies the effect of varying sparsity patterns | Moderate dimension, variable support | $n = 300$, $p = 1000$, support size varies, SNR grid, tFDR = 0.1, MC = 200 |
| **04** | MC: L-loop Strategies | Compares different L-calibration strategies | Moderate dimension | $n = 300$, $p = 1000$, $s = 10$, random support, SNR $ \in \{0.1, 1, 5\} $, tFDR = 0.1, MC = 200 |
| **05** | Memory-Mapped Single Run | Demonstrates memory-mapped workflows | Larger data and reduced RAM usage | Two scenarios: in-memory + solver mmap, and full mmap pipeline |
| **06** | MC: Memory-Mapped | Studies FDR/TPR under memory-mapped execution | Large-scale setting | $n = 5000$, $p = 10000$, SNR grid, mmap storage on disk |
| **07** | MC: Scalability | Reserved for future scalability experiments | Placeholder | To be completed |

---

## Folder contents

```txt
trex/
  ├── README.md
  ├── CMakeLists.txt
  ├── trex_sim_utils.hpp
  │
  ├── demo_trex_01_single_run/
  │   ├── demo_trex_01_single_run.cpp
  │   ├── README.md
  │   └── simulation_results/
  │       ├── d01_trex_basic_n150_p300.txt
  │       └── d01_trex_basic_n5000_p1000.txt
  │
  ├── demo_trex_02_mc_sim_fixed_support/
  │   ├── demo_trex_02_mc_sim_fixed_support.cpp
  │   ├── README.md
  │   └── simulation_results/
  │       ├── demo_trex_02_mc_sim_fixed_support_trex_results_n300_p1000_stagnation_window_7.txt
  │       └── demo_trex_02_mc_sim_fixed_support_trex_results_n300_p1000_stagnation_window_7.csv
  │
  ├── demo_trex_03_mc_sim_variable_support/
  │   ├── demo_trex_03_mc_sim_variable_support.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_04_mc_sim_lloop_strategies/
  │   ├── demo_trex_04_mc_sim_lloop_strategies.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_05_mmap/
  │   ├── demo_trex_05_mmap.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  ├── demo_trex_06_mc_sim_mmap/
  │   ├── demo_trex_06_mc_sim_mmap.cpp
  │   ├── README.md
  │   └── simulation_results/
  │
  └── demo_trex_07_mc_sim_scalablity/
      ├── demo_trex_07_mc_sim_scalablity.cpp
      ├── README.md
      └── simulation_results/
```

---

## Example results

The following example comes from **Demo 02**, which evaluates performance across several SNR levels for multiple base
solvers within the T-Rex framework. We depict here an exemplary result using the TLARS solver.

| SNR | FDR (TLARS) | TPR (TLARS) |
|-----|-------------|-------------|
| 0.1 | 0.0050      | 0.0030      |
| 0.5 | 0.0481      | 0.2600      |
| 1.0 | 0.0427      | 0.7360      |
| 2.0 | 0.0398      | 0.9845      |
| 5.0 | 0.0279      | 1.0000      |

### Interpretation

- FDR remains below the target level $\text{tFDR} = 0.1$ across all listed SNR values.
- TPR increases strongly as the signal becomes easier to detect.
- At very low SNR, recovery is difficult.
- At moderate to high SNR, the procedure recovers most or all active variables.

This is the expected pattern for a method that prioritizes controlled false discoveries while gaining power as signal strength increases.

---

## Building the demos

From the C++ workspace root:

```bash
cd TRexSelector_Simulations/cpp
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_PREFIX_PATH="<path-to-TRexSelector>/cpp/install"
cmake --build build/debug
```

The executables are written to:

```txt
build/debug/bin/
```

---

## Running a demo

For example, to run Demo 02:

```bash
./build/debug/bin/demo_trex_02_mc_sim_fixed_support
```

The program writes its output to the console and also stores result files in the corresponding demo folder.
From here the results can be used in a graphics library of choice, such as Python's `matplotlib`, R's `ggplot2`,
Cpp's `matplot++` or any other plotting tool.

Example output location:

```txt
demo_trex_02_mc_sim_fixed_support/simulation_results/
  ├── demo_trex_02_mc_sim_fixed_support_trex_results_n300_p1000_stagnation_window_7.txt
  └── demo_trex_02_mc_sim_fixed_support_trex_results_n300_p1000_stagnation_window_7.csv
```

### Output files

Depending on the demo, the `simulation_results/` folder may contain:

- **`.txt` files** for human-readable summaries,
- **`.csv` files** for plotting, post-processing, or comparison across methods.

---

## Notes for new users

- Start with **Demo 01** if you only want to verify that the pipeline works.
- Use **Demo 02** and **Demo 03** to study statistical behavior across repeated simulations.
- Use **Demo 05** and **Demo 06** if you are interested in larger problem sizes or disk-backed workflows.
- Check the local `README.md` inside each demo folder for scenario-specific details.

---

## References

- The main T-Rex Selector implementation and methodological background are documented in the parent project.
- See [TRexSelector/README.md](../../../README.md) for architecture, dependencies, and package-level information.

---

**Last updated**: 2026-07-01
