# Holographic MIMO Channel Estimation

Code accompanying the code-aided channel estimation framework for metasurface-based
Holographic MIMO (HMIMO) systems, covering both **Stacked Intelligent Metasurface (SIM)**
and **Beyond-Diagonal SIM (BD-SIM)** architectures.

## Papers

If you use this code, please cite the following works:

**Main reference (this repository):**

> R. C. G. Porto and R. C. de Lamare, "Code-Aided Channel Estimation for
> Metasurface-Based Holographic MIMO Systems," accepted at *IEEE VTC 2026 Spring*.
> Preprint: https://arxiv.org/abs/2606.22465

**Original iterative framework this work builds upon:**

> R. C. G. Porto and R. C. de Lamare, "Iterative Joint Channel Estimation and
> Detection for Coded Multi-RIS-Assisted Multi-Antenna Systems," *IEEE Transactions
> on Communications*, 2025. https://ieeexplore.ieee.org/abstract/document/11313562/

## ⚠️ Important note on this release

This code underwent a **substantial refactoring** after the results in the paper
were obtained. The original implementation used separate, duplicated scripts per
architecture and per number of metasurface layers (e.g., one file for SIM with 1
layer, another for SIM with 3 layers, another for BD with 5 layers, etc.). This
release consolidates all of that into a small set of unified, parameterized
functions (see `lib/`), dispatching on architecture (`"SIM"`/`"BD"`) and number of
layers (`R`) as configuration, rather than as separate code paths.

**Every refactored function was validated against the original implementation**
using two complementary approaches:

- **Bit-exact sanity checks**, where feasible, confirming identical outputs under
  matched random seeds.
- **Statistical equivalence tests** (mean/variance comparison and two-sample
  t-tests across N independent Monte Carlo trials), used whenever the refactoring
  changed the internal order of random number consumption. This was necessary in
  a few cases (e.g., the unified BD layer sampling), where the algorithm is
  provably identical but produces a different — statistically equivalent —
  realization sequence.

All validation scripts are included in `tests/`. That said, this was a large
reorganization done in a short time frame: **some edge cases or configurations may
still need minor adjustments**. If you run into unexpected behavior, please open
an issue — feedback is very welcome and will help improve the release for future
users.

## Repository structure

```
.
├── configs/          # System initialization (parameters, LDPC, geometry, etc.)
├── deps/              # Shared utilities (LDPC codec, channel generators, misc.)
├── lib/               # Core unified methods:
│   │                  #   - ce_lmmse.m   : channel estimator (SIM/BD, any R)
│   │                  #   - idd_mmse.m   : IDD detector (SIM/BD, any R)
│   │                  #   - alt_opt.m    : metasurface phase/unitary design (AO)
│   │                  #   - compute_wk.m : MMSE-SIC filter design
│   │                  #   - run_pipeline.m : full coarse→refined→perfect pipeline
├── experiments/       # Entry-point scripts reproducing the paper's results
├── tests/             # Statistical/sanity validation of the unified functions
├── figures/           # Generated figures (empty by default)
├── results/           # Generated .mat results (empty by default)
└── hmimo_channel_estimation.prj   # MATLAB project file
```

## Requirements

- MATLAB (developed and tested on R2024a or newer; earlier versions may work but
  are untested).
- Communications Toolbox (LDPC coding/decoding, QAM modulation).
- Signal Processing Toolbox.
- Parallel Computing Toolbox (optional — `parfor` is used in the experiment
  scripts for speed, but falls back to a regular loop if unavailable).

## How to use

1. Clone this repository.
2. Open `hmimo_channel_estimation.prj` in MATLAB — this sets up the correct path
   automatically.
3. Run one of the scripts in `experiments/`, for example:

   ```matlab
   hmimo_template
   ```

   This runs the full coarse → refined → perfect-CSI pipeline over a sweep of
   transmit power values, for a chosen architecture (`"SIM"` or `"BD"`) and number
   of metasurface layers `R`, and plots BER, sum-rate, and NMSE.

4. To change the scenario, edit the top of the experiment script:

   ```matlab
   Pw_dBm_vec      = 3:1:5;     % transmit power sweep (dBm)
   NumberOfPackets = [10 20 30];% Monte Carlo packets per power value
   n_layers        = 3;         % number of metasurface layers (R)
   arch            = "SIM";     % "SIM" or "BD"
   ```

5. All system parameters (antenna counts, LDPC code, geometry, noise figure,
   etc.) are set in `configs/init_system.m`.

### Running the validation tests

Before trusting results in a new environment, you can re-run the equivalence
tests in `tests/` to confirm the refactored functions behave as expected. These
compare the unified functions against the (now superseded) per-architecture,
per-layer implementations across many independent trials.

## Scope of this release

This repository currently contains the **HMIMO channel estimation code** (SIM and
BD-SIM architectures) corresponding to the VTC 2026 Spring paper above. 

## Contact

Roberto C. G. Porto
camara@ime.eb.br / camara2k@gmail.com

## License

MIT License.
