# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RAMSES-CPP is a modern C++17 port of the legacy [RAMSES-2025](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Fortran Adaptive Mesh Refinement (AMR) astrophysics code. The primary goal is **strict binary parity** with the original Fortran snapshots while offering a modular, modern C++ architecture. The legacy Fortran source lives in `./legacy/` and is the reference for porting logic.

## Build System

CMake >= 3.15 is required. The build produces three executables (`ramses_1d`, `ramses_2d`, `ramses_3d`) and a verification tool (`verify_ref`).

```bash
mkdir -p build && cd build
cmake .. -DRAMSES_NDIM=3 -DRAMSES_USE_MPI=OFF -DRAMSES_USE_MHD=ON -DRAMSES_USE_RT=ON
make -j$(nproc)
```

Key CMake flags:
- `RAMSES_NDIM` — Dimensionality: `1`, `2`, or `3` (default `3`)
- `RAMSES_USE_MHD` — Enable magnetohydrodynamics (`ON`/`OFF`)
- `RAMSES_USE_RT` — Enable radiative transfer (`ON`/`OFF`)
- `RAMSES_USE_MPI` — Enable MPI (`ON`/`OFF`, default `OFF`)
- `RAMSES_PRECISION` — Float precision: `4` or `8` (default `8`)
- `RAMSES_NENER` — Number of non-thermal energy groups
- `RAMSES_NPSCAL` — Number of passive scalars
- `RAMSES_NGROUPS` / `RAMSES_NIONS` — RT photon groups / ion species
- `RAMSES_LONGINT` — Use 64-bit grid IDs (`ON`/`OFF`)

These flags become compile-time `#define` macros (e.g., `NDIM=3`, `MHD`, `RT`). The `NDIM` macro is defined per target so all three dimensionality libraries are built simultaneously from the same sources.

## Running Simulations

Simulations take a RAMSES `.nml` namelist file as the sole argument:

```bash
./build/ramses_3d namelist/sedov3d.nml
# MPI
mpirun -np 8 ./build/ramses_3d namelist/sedov3d.nml
```

Pre-configured namelists live in `namelist/`. Each test case in `tests/` has its own `.nml` alongside a `config.txt` that records the `FLAGS` used for compilation (e.g., `FLAGS: NDIM=1 SOLVER=hydro`).

## Testing

**Always use `timeout` — simulations can stall indefinitely due to numerical instability.**

```bash
export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu

# Run all tests
cd tests && timeout 10m ./run_test_suite.sh

# Run a specific test category
cd tests && timeout 10m ./run_test_suite.sh -t hydro
cd tests && timeout 10m ./run_test_suite.sh -t mhd

# Run specific test numbers
cd tests && timeout 10m ./run_test_suite.sh -t 3-5,10

# Parallel run (4 MPI ranks)
cd tests && timeout 10m ./run_test_suite.sh -p 4
```

The test runner reads each test's `config.txt`, re-invokes `cmake` with the appropriate flags, runs the simulation, and calls the `plot-*.py` script which prints `PASSED` or `FAILED`. Results are written to `tests/test_suite.log` and a `test_results.pdf`.

**Test file constraints:** Do NOT edit `*.nml`, `*-ref.dat`, or `plot-*.py` files inside `tests/`. Shell scripts (`*.sh`) in `tests/` may be modified.

### Snapshot Parity Verification

```bash
cd build && ./verify_ref <reference_snapshot> <local_snapshot>
# e.g.:
./verify_ref tests/hydro/sod-tube/output_00001/hydro_00001.out00001 \
             /tmp/output_00001/hydro_00001.out00001
```

## Architecture

### Core Data Structure: `AmrGrid`

`AmrGrid` (`include/ramses/AmrGrid.hpp`) is the central data store for the entire AMR octree. It owns all grid and particle arrays as flat `std::vector`s with Fortran-style 1-based indexing. Key conventions:

- Coarse cells (level 0) occupy indices `1..ncoarse`; refined cells start at `ncoarse+1`
- `headl(icpu, ilevel)` / `taill` / `numbl` — linked lists of grids per (rank, level)
- `son[igrid]` — index of the first child cell; `father[igrid]` — parent cell index
- `nbor[igrid * 2*NDIM + iface]` — neighbor grid indices for ghost zone exchange
- `xg[igrid * 3]` — oct center coordinates (3D storage even for 1D/2D; unused dims default to `0.5 * boxlen`)
- `uold_vec` / `unew_vec` — conservative hydro variables, layout: `[nvar][ncell]`
- Particles: `xp`, `vp`, `mp`, `tp`, `zp`, `idp`, `family`, `tag` — families are `FAM_DM=1`, `FAM_STAR=2`, `FAM_SINK=3`, `FAM_TRACER=6` (see `Types.hpp`)

### `Simulation` and the Time-Step Loop

`Simulation` (`src/Simulation.cpp`) drives everything:
1. `initialize(nml_path)` — parses the `.nml` file via `Config`, instantiates all solvers through `SolverFactory`, and builds the initial AMR tree via `Initializer` inside the refinement loop.
2. `run()` — top-level time loop calling `amr_step(0, dt)`.
3. `amr_step(ilevel, dt)` — recursive sub-cycling: applies hydro/MHD/RT/gravity at each level before recursing to `ilevel+1`. Calls `TreeUpdater` to refine/coarsen after each level step.

### Solver Pattern

All physics modules are polymorphic classes instantiated by free functions in `SolverFactory` (`src/SolverFactory.cpp`). Each solver takes `AmrGrid&` and `Config&` in its constructor and is owned by `Simulation` via `std::unique_ptr`. To add a new physics module: implement a class derived from the relevant base, register a factory function in `SolverFactory.hpp/.cpp`, and add the `std::unique_ptr` member to `Simulation`.

Key solvers:
- `HydroSolver` — MUSCL-Hancock Godunov scheme; slope limiting via `SlopeLimiter`; Riemann solvers in `RiemannSolver`
- `MhdSolver` — Constrained-transport MHD on staggered grid (compile with `-DRAMSES_USE_MHD=ON`)
- `RhdSolver` — Relativistic hydrodynamics with Newton-Raphson primitive recovery and TM EOS
- `RtSolver` / `RtChemistry` — M1 radiative transfer (compile with `-DRAMSES_USE_RT=ON`)
- `PoissonSolver` — Multi-grid self-gravity with comoving support
- `ParticleSolver` — N-body + tracer particles with trilinear interpolation
- `StarSolver` / `FeedbackSolver` — Star formation (Poisson statistics) and SN feedback
- `SinkSolver` — Sink particles for sub-grid accretion
- `TurbulenceSolver` — Stochastic spectral driving in Fourier space
- `ClumpFinder` — On-the-fly density peak / saddle-point structure identification

### I/O and Parity

`RamsesWriter` / `RamsesReader` (`src/`) produce binary snapshots in **Fortran unformatted record** format using 4-byte record markers, exactly matching legacy RAMSES output. The snapshot level indexing uses `il + 1` (Level 0 in C++ → Level 1 in file). `Hdf5Writer` provides a parallel HDF5 alternative. Visualization uses `tests/visu/visu_ramses.py` (set `PYTHONPATH` as above).

### Configuration Parsing

`Config` (`src/Config.cpp`) parses RAMSES `.nml` files. It supports standard `&BLOCK / key=value /` syntax, Fortran array keys (`key(1:2)=val1,val2`), and case-insensitive lookups. Access values via `config.get_int("RUN_PARAMS", "nstepmax", 10)`, etc.

### MPI

`MpiManager` is a singleton. `LoadBalancer` distributes octs using Hilbert-curve keys (`Hilbert.cpp`). Build with `-DRAMSES_USE_MPI=ON` to activate; otherwise all MPI calls are no-ops.

### Custom Patches

C++ source files placed in `patch/` are automatically discovered by CMake and compiled in with `-DRAMSES_USE_CUSTOM_PATCHES`. This is the mechanism for porting legacy Fortran override files (see `tools/ramses-patch-porter`).

## Code Style & Porting Patterns

- Translate Fortran idioms into safe, modern C++ patterns. Avoid C-style casts, unmanaged raw arrays (prefer `std::vector`), and magic numbers.
- Do not bypass the type system or disable warnings. Use explicit instantiation and strong typing.
- When porting logic, always cross-reference the original Fortran in `./legacy/`.

## Mandatory Documentation Updates

After any significant code change, update:
- `PORTING_HISTORY.md` — log the phase, what changed, and why
- `docs/` — update the relevant architecture/usage/physics markdown guide
- `README.md` — if build instructions or project status changed

After updating docs, you are **explicitly expected to `git commit` and `git push`** the changes.
