# T-Rex+DA: Dependency-Aware T-Rex — Python Demos

Python ports of the Dependency-Aware (DA) T-Rex demos. They mirror the R suite
in `R/trex_selector_methods/trex_da/` and the C++ suite in
`cpp/trex_selector_methods/trex_da/`.

Shared infrastructure:

- [dgp_generators.py](dgp_generators.py) — the DA data-generating processes
  (AR(1); equicorrelated/BT factor model; NN/MA(kappa); AR(1)-block, with and
  without a white block; heavy-tailed block, with and without white; prior-groups),
  all SNR-calibrated with 0-based support.
- [da_sim_common.py](da_sim_common.py) — support generators, a DA selector
  builder, a generic parallel Monte Carlo runner, and table/matrix printers.

The folder layout mirrors the C++ suite: one subfolder per demo (each with its
own `README.md` and `simulation_results/`), plus the two suite-level shared
modules above.

```txt
trex_da/
  ├── README.md
  ├── dgp_generators.py       ┐  suite-level shared modules
  ├── da_sim_common.py        ┘  (imported by every demo)
  ├── demo_trex_da_01_ar1/
  │   ├── demo_trex_da_01_ar1.py
  │   ├── README.md
  │   └── simulation_results/
  ├── demo_trex_da_02_equi/
  │   └── ...
  └── demo_trex_da_08_groups/
```

## Demos

| Folder | Description | cpp counterpart |
|---|---|---|
| [demo_trex_da_01_ar1](demo_trex_da_01_ar1/) | AR(1) Toeplitz: single run; SNR sweep (Random vs CappedSpread); rho sweep; 2D gap x rho (kappa boundary) | `demo_trex_da_01_mc_sim_ar1` |
| [demo_trex_da_02_equi](demo_trex_da_02_equi/) | Equicorrelated DGP (BT factor model, `rho_within == rho_between`), EQUI method | `demo_trex_da_02_mc_sim_equi_and_bt` |
| [demo_trex_da_03_nn](demo_trex_da_03_nn/) | Nearest-Neighbours / MA(kappa) banded DGP: single run; SNR/rho/kappa sweeps; 2D kappa x rho and SNR x rho | `demo_trex_da_03_mc_sim_nn` |
| [demo_trex_da_03b_nn_ar](demo_trex_da_03b_nn_ar/) | NN method on AR(1) data (method-mismatch): single run; SNR/rho sweeps; 2D SNR x rho | `demo_trex_da_03b_mc_sim_nn_ar` |
| [demo_trex_da_04_bt_ar1_block](demo_trex_da_04_bt_ar1_block/) | AR(1)-block DGP (p = M*Q): SNR/rho/Q/M sweeps, each with a linkage loop | `demo_trex_da_04_mc_sim_bt_ar1_block` |
| [demo_trex_da_05_bt_ar1_block_white](demo_trex_da_05_bt_ar1_block_white/) | AR(1)-block + white-noise (fixed p_total): SNR/rho/Q/M sweeps with linkage loop | `demo_trex_da_05_mc_sim_bt_ar1_block_sweeps` |
| [demo_trex_da_06_bt_ht_block](demo_trex_da_06_bt_ht_block/) | Heavy-tailed Student-t(3) blocks, Gaussian vs heavy noise: SNR/rho/Q/M/tFDR/linkage sweeps | `demo_trex_da_06_mc_sim_bt_ht_block_sweeps` |
| [demo_trex_da_07_bt_ht_block_white](demo_trex_da_07_bt_ht_block_white/) | Heavy-tailed blocks + heavy-tailed white, Gaussian vs heavy noise: SNR/rho/Q/M/tFDR sweeps | `demo_trex_da_07_mc_sim_bt_ht_block_white` |
| [demo_trex_da_08_groups](demo_trex_da_08_groups/) | Prior-groups method (`DAMethod.PRIOR_GROUPS`): 3-level hierarchy, SNR sweep | `demo_trex_da_08_mc_sim_groups` |

## Running

```bash
python Python/trex_selector_methods/trex_da/demo_trex_da_01_ar1/demo_trex_da_01_ar1.py
```

Each demo exposes `RUN_PART_*` flags at the top; a light default part runs
out of the box and the heavy Monte Carlo sweeps are opt-in (num_MC=200 each).
Reduce `num_MC` to preview. Each demo inserts both its own directory and its
parent suite directory onto `sys.path`, so the shared `dgp_generators` /
`da_sim_common` modules resolve from the nested location and are re-importable
by the `spawn`-launched workers — so demos must be run as files, not piped via
`python -`. Each demo folder has its own `README.md`.

Indices are **0-based** (Python convention); the R counterparts are 1-based.

**Last updated**: 2026-07-08
