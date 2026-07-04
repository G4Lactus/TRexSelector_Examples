# Demo 05: Memory-Mapped Data — Single Run

## Purpose

Demonstrate **memory-mapped (mmap) I/O workflows** for large datasets that do not fit in RAM. This demo shows how the T-Rex selector can handle design matrices $\mathbf{X}$ and dummy matrices $\mathbf{D}$ via disk-based storage with minimal performance overhead.

---

## Data Generation Parameters

- **Sample size**: $n = 5000$ (large)
- **Number of features**: $p = 10000$ (very high-dimensional)
- **True support**: $\mathcal{S}^* = \{$first 20 features$\}$
- **True coefficients**: Uniform $\pm 1$
- **SNR**: $1.0$ or $2.0$
- **DGP**: $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\epsilon}$, $\boldsymbol{\epsilon} \sim N(0, I_n)$

**Data size estimate**: 
- $\mathbf{X}$: $5000 \times 10000 \times 8$ bytes $\approx 400$ GB (impractical to hold in RAM)
- Mmap allows working with disk-resident files as if they were in-memory arrays

---

## Scenarios

### Scenario A: In-Memory $\mathbf{X}$, Mmap $\mathbf{D}$ (Dummies)
- Design matrix $\mathbf{X}$ stored in RAM
- Dummy matrices $\mathbf{D}$ created and stored on disk (mmap'd)
- Solvers read dummy columns on-demand via page fault mechanism

**Use case**: Medium-scale problems where $\mathbf{X}$ fits but $\mathbf{D}$ does not

### Scenario B: Fully Mmap'd ($\mathbf{X}$ and $\mathbf{D}$)
- Both $\mathbf{X}$ and $\mathbf{D}$ stored on disk
- All matrix access goes through mmap interface
- Maximum memory efficiency

**Use case**: Very large-scale problems; worst-case IO overhead

---

## Control Parameters

```
K = 20                           # Random experiments per T-loop iteration
max_dummy_multiplier = 10        # Max L = 10p
use_max_T_stop = true            # Cap T ≤ ceil(n/2)
tloop_stagnation_stop = true     # Early exit
opt_threshold = 0.75             # Optimization grid point
target_fdr = 0.1                 # Target FDR
use_memory_mapping = VARIES      # Main variable: toggle mmap on/off
mmap_data_dir = "/tmp/trex_mmaps/" # Disk location for mmap files
```

---

## Output Files

### Scenario A Output
**`demo_trex_05_mmap_scenario_A_n5000_p10000_snr1p0.txt`**

Contains selection results and timing:
```
=== T-Rex Selector: Mmap Scenario A ===
Design matrix (X) location: RAM
Dummy matrices (D) location: Disk (/tmp/trex_mmaps/)

Selected variables: ...
True positives: ...
FDR: ...

Timing:
  Data generation: 12.34 s
  T-Rex execution (mmap D): 45.67 s
  Total: 57.98 s
```

### Scenario B Output
**`demo_trex_05_mmap_scenario_B_n5000_p10000_snr2p0.txt`**

Similar format; timing will be higher due to $\mathbf{X}$ page faults.

---

## Running the Demo

```bash
./build/debug/bin/demo_trex_05_mmap
```

**Computation time**: 5–15 minutes depending on disk speed and mmap page fault rate.

**Disk space**: The demo creates temporary mmap files in `/tmp/trex_mmaps/` (or configured directory). Ensure ~1–2 GB free space.

---

## Expected Results

- **Selection quality**: Identical FDR/TPR regardless of mmap setting (FDP estimation is unaffected by storage medium)
- **Timing overhead**: Scenario B typically runs 1.5–3× slower than Scenario A due to $\mathbf{X}$ page faults
- **Memory footprint**: RAM usage is minimal when mmap is fully enabled (Scenario B)

---

## Interpretation Guide

**What to look for:**
1. **Correctness**: Does mmap produce the same selection as in-memory data?
   - Expected: Yes, identical results
2. **Performance trade-off**: How much faster is Scenario A vs. B?
   - Expected: A is 1.5–3× faster
3. **Disk IO statistics**: How many page faults occur?
   - Expected: Thousands to tens of thousands (visible in system monitoring)

**Practical guidance:**
- Use mmap when RAM is insufficient
- For optimal performance: keep $\mathbf{X}$ in RAM if possible; only mmap $\mathbf{D}$
- Disk speed significantly impacts runtime (SSD >> HDD)

---

## Cleanup

After the demo, remove temporary mmap files:
```bash
rm -rf /tmp/trex_mmaps/
```

---

**Last updated**: 2026-07-01
