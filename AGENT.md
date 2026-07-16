# AGENT.md

This file provides guidance to AI coding agents when working with code in this repository.

## Project Overview

RAMSES-CPP is a modern C++17 port of the legacy [RAMSES-2025](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Fortran Adaptive Mesh Refinement (AMR) astrophysics code. The primary goal is **strict binary parity** with the original Fortran snapshots while offering a modular, modern C++ architecture. The legacy Fortran source lives in `./legacy/` and is the reference for porting logic.

## Constraints

- **Do NOT edit any test files** (`*.nml`, `*-ref.dat`, `plot-*.py` inside `tests/`). Only shell scripts (`*.sh`) in `tests/` may be modified.

## Build System

CMake >= 3.15 is required. The build produces three optimized executables (`ramses_1d`, `ramses_2d`, `ramses_3d`).

```bash
mkdir -p build && cd build
# Release builds are 9.2x faster than Debug (default is Release)
cmake .. -DNDIM=3 -DMPI=OFF -DSOLVER=mhd -DRT=ON
make -j$(nproc)
```

Key CMake flags:
- `NDIM` — Dimensionality: `1`, `2`, or `3` (default `3`)
- `SOLVER` — Solver type: `hydro`, `mhd`, or `rhd` (default `hydro`)
- `RT` — Enable radiative transfer (`ON`/`OFF`)
- `MPI` — Enable MPI (`ON`/`OFF`, default `OFF`)
- `NPRE` — Float precision: `4` or `8` (default `8`)
- `NENER` — Number of non-thermal energy groups
- `NPSCAL` — Number of passive scalars
- `NGROUPS` / `NIONS` — RT photon groups / ion species
- `LONGINT` — Use 64-bit grid IDs (`ON`/`OFF`)

These flags match legacy Makefile options. The build system only compiles the target dimensionality specified by `NDIM` to optimize build times.

## Running Simulations

Simulations take a RAMSES `.nml` namelist file as the sole argument. **ALWAYS use `timeout`** to prevent indefinite stalls due to numerical instability:

```bash
# Set PYTHONPATH for visualization/analysis scripts
export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu

./build/ramses_3d namelist/sedov3d.nml

# MPI
mpirun -np 8 ./build/ramses_3d namelist/sedov3d.nml
```

Pre-configured namelists live in `namelist/`. Each test case in `tests/` has its own `.nml` alongside a `config.txt` that records the `FLAGS` used for compilation.

## Testing

```bash
# Run all tests (with 10m timeout)
cd tests && timeout 10m ./run_test_suite.sh

# Run a specific test category
cd tests && timeout 10m ./run_test_suite.sh -t hydro
cd tests && timeout 10m ./run_test_suite.sh -t mhd

# Run a single test specifically (e.g. advect1d is test 1) and keep/save results data
cd tests && ./run_test_suite.sh -t 1 -k
```

### Test Suite Options
* `-t <category>`: Run a specific test category (e.g. `hydro`, `mhd`) or test range/indices (e.g., `-t 3-5,10` or `-t 1` for `advect1d`).
* `-p <N>`: Run in parallel using MPI on `<N>` processes.
* `-r`: Run tests with restart verification.
* `-c`: Clean build and delete intermediate build files before starting.
* `-k`: Keep/save results data (do not delete output folders like `output_*` after comparison finishes).
* `-s`: Run the test suite with coverage collection (`gcov`).

### Snapshot Parity Verification
The primary goal is **binary parity** with legacy RAMSES. Output snapshots are automatically compared to legacy references by the Python scripts in `tests/visu/` when running the test suite.


## Architecture

### Core Data Structure: `AmrGrid`

`AmrGrid` (`include/ramses/AmrGrid.hpp`) is the central data store for the entire AMR octree.
- **Dynamic Resizing:** Grid storage grows automatically via `resize_grids()` to prevent overflow during refinement.
- **1-Based Indexing:** Uses Fortran-style 1-based indexing for parity. Coarse cells are Level 0; refined grids are Level 1+.
- **Linked Lists:** `headl(icpu, ilevel)` tracks grids; `son` and `father` arrays define the hierarchy.

### Solver Pattern & Phase 45 Updates

All physics modules are polymorphic classes instantiated by `SolverFactory`. 
**Phase 45** integrated missing Riemann solvers:
- `RiemannSolver::solve_acoustic` — Exact acoustic solver.
- `RiemannSolver::solve_godunov_nr` — Exact multi-dimensional Godunov solver (Newton-Raphson).
- **Slope Alignment:** `slope_type = 0` now correctly maps to zero slopes (first-order).

Key solvers:
- `HydroSolver` — MUSCL-Hancock Godunov scheme; supports `exact`, `hllc`, `llf`, and `acoustic` Riemann solvers.
- `MhdSolver` — Constrained-transport MHD on staggered grid.
- `ParticleSolver` — N-body + tracer particles (FAM_TRACER=6).

## Code Style & Porting Patterns

- **Heritage:** Always cross-reference logic with the legacy Fortran in `./legacy/`.
- **Strict Parity:** Do not change numerical algorithms or output formats unless specifically fixing a divergence from legacy.
- **Documentation:** After significant changes, update `PORTING_HISTORY.md` and relevant docs in `docs/`.

## Mandatory Updates

After significant work, update:
- `PORTING_HISTORY.md` — log the phase and changes.
- `README.md` — if build/status changed.
- `AGENT.md` (this file) — if conventions or architecture shifted.

You are **explicitly expected to `git commit` and `git push`** updates to documentation.
