---
layout: default
title: Usage
---

# Usage Guide

RAMSES-CPP is fully compatible with the legacy RAMSES ecosystem, including namelists and visualization tools.

## Running Simulations

To start a simulation, choose the executable matching your problem dimensionality and provide a standard RAMSES namelist file.

```bash
# Standard 3D execution
./ramses_3d path/to/namelist.nml

# MPI-scaled execution (e.g., 8 ranks)
mpirun -np 8 ./ramses_3d path/to/namelist.nml
```

### Namelist Compatibility
The project uses the standard RAMSES namelist format. Key blocks include:
- `&RUN_PARAMS`: Control the simulation lifecycle (`hydro`, `poisson`, `ncontrol`, etc.).
- `&AMR_PARAMS`: Grid resolution and refinement (`nx,ny,nz`, `levelmin`, `levelmax`).
- `&HYDRO_PARAMS`: Fluid physics settings (`gamma`, `riemann`, `slope_type`).
- `&SF_PARAMS`: Star formation settings (`n_star`, `eps_star`, `m_star`).
- `&FEEDBACK_PARAMS`: Supernova and metal parameters (`yield`, `eta_sn`).
- `&CLUMP_PARAMS`: Structure identification thresholds.
- `&RT_PARAMS`: Radiation groups and parameters.

## Star Formation and Feedback

To enable star formation and stellar feedback:

```fortran
&RUN_PARAMS
hydro=true
poisson=true
star=true
/

&SF_PARAMS
n_star=0.1      ! SF density threshold [H/cc]
eps_star=0.01   ! Efficiency per free-fall time
m_star=1.5      ! Base star mass [Code units]
/

&FEEDBACK_PARAMS
eta_sn=0.1      ! SN energy fraction
yield=0.02      ! Metal yield
/
```

## Structure Identification (ClumpFinder)

Enable clump finding in the snapshot dump:

```fortran
&RUN_PARAMS
clumpfind=true
/

&CLUMP_PARAMS
density_threshold=100.0
/
```

## Visualization and Analysis

RAMSES-CPP maintains **Bit-Perfect Parity** in its binary outputs. This means you can use legacy `visu_ramses.py` or OSIRIS directly.

### Python Plotting
1. **Set PYTHONPATH:**
   ```bash
   export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu
   ```
2. **Execute Plotting Script:**
   ```bash
   python3 tests/hydro/advect1d/plot-advect1d.py
   ```

## Automated Test Suite

A comprehensive test suite is included to ensure physics and parity integrity.

### Mandatory Watchdog Rule
Always use the `timeout` command to prevent hangs during development:
```bash
cd tests
timeout 10m ./run_test_suite.sh -t hydro
```

## Snapshot Verification

Use the `verify_ref` tool to compare your results against legacy Fortran reference snapshots:
```bash
./verify_ref path/to/ref_amr_00001.out00001 path/to/local_amr_00001.out00001
```

---
*Pro Tip: Check the `namelist/` directory for pre-configured benchmark examples!* 🛰️✨
