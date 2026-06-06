---
layout: default
title: Usage
---

# Simulation Usage Guide

This document describes how to execute simulations in RAMSES-CPP, configure namelist parameters, and verify outputs.

---

## 🏃 Running Simulations

RAMSES-CPP generates dimension-specific executables named `ramses_1d`, `ramses_2d`, and `ramses_3d` based on the configured `NDIM` CMake flag.

### 1. Serial Execution
To run a serial simulation, pass a standard RAMSES `.nml` namelist file as the sole argument:
```bash
./build/ramses_3d namelist/sedov3d.nml
```

### 2. MPI Parallel Execution
To run on multiple processors, use the standard `mpirun` wrapper:
```bash
mpirun -np 8 ./build/ramses_3d namelist/sedov3d.nml
```

> [!WARNING]
> **Safety Watchdog:** To prevent simulations from hanging indefinitely due to physics instabilities or tiny timesteps during development, always run tests or debug jobs with the `timeout` command:
> `timeout 10m ./build/ramses_3d namelist/sedov3d.nml`

---

## ⚙️ Namelist Parameters

Simulation runs are configured using the standard Fortran namelist format, organized into parameter blocks. Key blocks include:

### `&RUN_PARAMS`
Controls the physical modules and solvers active in the run.
* `hydro=true`: Enables the hydrodynamics/MHD fluid solver.
* `poisson=true`: Enables the self-gravity Poisson solver.
* `particles=true`: Enables N-body dark matter/stellar particles.
* `pic=true`: Enables Particle-Mesh calculations.

### `&AMR_PARAMS`
Governs the grid resolution, refinement criteria, and domain size.
* `nx, ny, nz`: Grid cell counts at the root level (Level 0).
* `levelmin`: The base grid level (e.g., `4` corresponds to a $2^4 = 16$ grid).
* `levelmax`: The maximum allowed refinement level.
* `boxlen`: Size of the simulation domain (usually in code units).

### `&HYDRO_PARAMS`
Physical and numerical fluid dynamics configurations.
* `gamma`: Adiabatic index (e.g., `1.4` for diatomic gas, `1.667` for monoatomic).
* `courant_factor`: The CFL timestep restriction multiplier (default `0.8`).
* `riemann`: The Riemann solver to use (`hllc`, `llf`, `exact`, etc.).
* `slope_type`: Slope limiter for MUSCL reconstruction (`1`=MinMod, `2`=MC, `3`=Superbee).

---

## 📊 Visualizing Results

Since RAMSES-CPP preserves **Bit-Perfect Binary Parity** in its output structures, you can analyze and visualize output snapshots using existing legacy python scripts (such as the standard `visu_ramses.py` or OSIRIS pipelines).

### Python Setup Example
To plot grid slices using the visualization modules in the test suite:
1. Export the python path containing the visualization libraries:
   ```bash
   export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu
   ```
2. Execute the python plotting scripts:
   ```bash
   python3 tests/hydro/advect1d/plot-advect1d.py
   ```
This will parse the raw binary snapshot and generate a high-fidelity PDF showing the fluid variables alongside reference comparisons.
