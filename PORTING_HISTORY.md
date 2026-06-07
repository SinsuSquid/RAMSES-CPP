# 🕰️ RAMSES-CPP Porting History

This document tracks the milestones, architecture updates, and physics solver integration details during the migration from legacy RAMSES (Fortran) to the modern C++17 engine.

---

## 🚩 Phase 48: Advect1d AMR & Solver Realignment (In Progress) 🚀✨
**Commit:** `bce0823`
* **Initial Grid Refinement:** Extended the C++ further refinement iterations in [Simulation.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/Simulation.cpp#L219) from `lmax` to `lmax + 2`. This aligned C++ with the combined refinement sweeps of Fortran's `init_refine` and `init_refine_2`, achieving bit-perfect grid matching in Snapshot 1.
* **Grid Deletion Correction:** Corrected `remove_grid_fine` in [TreeUpdater.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/TreeUpdater.cpp#L140) to check if the father cell at level `ilevel - 1` is unflagged (`flag1 == 0`), matching Fortran's grid deletion logic.
* **Refinement Rules subcycling correction:** Aligned C++ subcycle checking in [TreeUpdater.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/TreeUpdater.cpp#L477) to evaluate the current level timing compared to the parent level subcycle factor.
* **Ultrabee Limiter formulation matching:** Removed the non-standard C++ safety fallback check in `HydroSolver::compute_slopes` for the Ultrabee limiter (`slope_type == 5`), which was overriding compressive slopes with minmod slopes. This brought Snapshot 2 grid size from 110 cells down to 104 cells, and aligned C++ time steps perfectly with Fortran.
* **Interpolation Modes & Limiters:** Implemented options for conserved variables (Mode 0), internal energy (Mode 1), and primitive variables (Mode 2) interpolation at coarse-fine boundaries, matching legacy `interpol_var` behavior, and added support for the combined central/MonCen limiter (`interpol_type = 4`).

---

## 🚩 Phase 47: CMake Parameter Alignment, verify_ref Cleanup, and Ninja Support (Completed) ✨🎯
**Commit:** `1277719`
* **CMake Variable Alignment:** Mapped and standardized all C++ build-system cache flags in [CMakeLists.txt](file:///home/bgkang/Projects/RAMSES-CPP/CMakeLists.txt) to match legacy Fortran variables exactly.
* **verify_ref Executable Cleanup:** Removed the obsolete `verify_ref` executable target and deleted `src/verify_ref.cpp`.
* **MIT License:** Integrated the MIT License under the copyright of `SinsuSquid`.
* **Ninja Build Support:** Updated [run_test_suite.sh](file:///home/bgkang/Projects/RAMSES-CPP/tests/run_test_suite.sh) to detect and use the `Ninja` generator automatically if available, yielding faster compile/re-compile times during local regression testing.
* **Advection Test Parity Diagnosis:** Successfully diagnosed the root cause of the advection test snapshot mismatch, identifying that the base-level timestep collapses to `1e30` when all cells at `levelmin` are refined, satisfying the premature dump condition.

---

## 🚩 Phase 46: Recursive Subcycling & Strict Initialization Alignment (Completed) ✨🎯
**Commit:** `a7431e3` & `5e00381`
* **Recursive Time Integration:** Refactored the time-stepping loop in [Simulation::amr_step](file:///home/bgkang/Projects/RAMSES-CPP/src/Simulation.cpp#L126) to mirror legacy `amr_step.f90` exactly. Timestep calculation (`newdt_fine`) and flagging are now integrated recursively, with levels managing their own timing via `dtnew_` and `dtold_` arrays.
* **Double-Fluxing Fix:** Implemented a robust `is_refined` check in [HydroSolver.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp) combining grid-based and cell-level flags to prevent double-counting fluxes at coarse-fine boundaries.
* **make_grid_fine Neighbor Lookup:** Corrected neighbor lookup in [TreeUpdater.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/TreeUpdater.cpp) to use proper coordinate math for coarse cells and sibling/cousin indexing for fine cells, preventing garbage memory interpolation.
* **Localized Initialization:** Implemented level-by-level initialization matching legacy `init_refine.f90`, ensuring AMR grids only grow where gradients exist.
* **Refluxing Parity:** Restored the $1 / 2^{\text{NDIM}}$ scaling factor for flux correction.

---

## 🚩 Phase 44: Slope Limiters, HLLC Energy, and Diagnostics (Completed) ✨
**Commit:** `778579f`
* **New Slope Limiters:** Implemented legacy slope types 6 (density-only centered), 7 (van Leer harmonic mean), and 8 (van Leer generalized MonCen) in [HydroSolver::compute_slopes](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp).
* **HLLC Star-Region Energy:** Replaced approximate star state calculations in [RiemannSolver::solve_hllc](file:///home/bgkang/Projects/RAMSES-CPP/src/RiemannSolver.cpp) with Toro's exact Rankine-Hugoniot jump condition formulations.
* **Grid Collapse Guard:** Added `NDIM` dimensionality guards to `ensure_ref_rules` to prevent AMR grid collapse in 1D and 2D.
* **Memory Diagnostics:** Restricted density checks to active coarse and valid parent-level grids, eliminating spurious `min_rho = 0` prints.

---

## 🚩 Phase 43: Advect1d AMR Parity — Refluxing & Refinement Rules (Completed) ✨
* **Refinement Rules (`ensure_ref_rules` & `authorize_fine`):** Ported nested grid verification logic from legacy `flag_utils.f90` and `virtual_boundaries.f90` to C++ [TreeUpdater](file:///home/bgkang/Projects/RAMSES-CPP/src/TreeUpdater.cpp). This checks that every refined cell at `ilevel` is surrounded by valid parent neighbors, zeroing violating refinement flags to enforce the strict 1-level difference rule.
* **`cell_levels` Cache Correction:** Adjusted the `cell_levels` cache in [HydroSolver.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp) to iterate through `ncpu + nboundary` (instead of just `ncpu`), resolving stride/offset issues when boundary domains are configured.
* **Refluxing Factor Alignment:** Aligned coarse-fine flux corrections (refluxing) in `HydroSolver.cpp` and `RhdSolver.cpp` with legacy physics by applying the correct $1 / 2^{\text{NDIM}}$ scaling factor.

---

## 🚩 Phase 42: AMR Refinement & Non-Thermal Energy (Completed) ✨
* **Refinement Flag Retainment (`test_flag`):** Updated [TreeUpdater::flag_fine](file:///home/bgkang/Projects/RAMSES-CPP/src/TreeUpdater.cpp) to check if a cell's children are already refined or flagged, preventing incorrect grid collapses.
* **Gradient-Based Flagging:** Enabled checks on density, pressure, velocity, and magnetic field gradients in `TreeUpdater::flag_fine`.
* **Non-Thermal Energy Initialization:** Corrected configuration parsing for multi-dimensional lists (e.g. `prad_region(1,1)`) in [Initializer.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/Initializer.cpp) and integrated non-thermal energy subtraction in [HydroSolver.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp) to properly compute internal energy.

---

## 🚩 Phase 41: Hydro Stability & Riemann Solver Alignment (Completed) ✨
* **Trace Step Flooring:** Implemented density/pressure flooring and slope clamping in [HydroSolver::trace](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp) to prevent numerical sound speed explosions.
* **Global CFL Timestep:** Corrected CFL calculation to scan all active refinement levels.
* **MUSCL-Hancock Predictor:** Integrated the half-step predictor calculation in `trace` for full second-order accuracy.
* **Exact Rankine-Hugoniot Solver:** Ported the iterative Newton-Raphson exact Riemann solver in [RiemannSolver::solve_godunov_nr](file:///home/bgkang/Projects/RAMSES-CPP/src/RiemannSolver.cpp) and aligned HLLC/LLF wave speed formulations with legacy Fortran.

---

## 🚩 Phase 40: MHD Module Stabilization & Output Compliance (Completed) ✨
**Commit:** `ee54e15` & `a7431e3`
* **MHD Crash Elimination:** Solved `std::out_of_range` exceptions in [MhdSolver::gather_stencil](file:///home/bgkang/Projects/RAMSES-CPP/src/MhdSolver.cpp) by replacing hardcoded 3D constant arrays with NDIM-aware stencil arrays.
* **Output Format Parity:** Standardized output files to always write all 3 velocity and 6 B-field components (setting unused components to 0.0) regardless of dimensionality, conforming to the legacy RAMSES file format.
* **Build System Optimization:** Re-architected CMake and `run_test_suite.sh` to compile only the targeted dimensionality (`RAMSES_NDIM`) and cached compilation flags to speed up test iteration times.

---

## 🚩 Phase 39: Solver Stability Investigation & Verification (Completed) ✨🎯
* **Passive Scalar Fixes:** Resolved Mach 2 mixing-scalar timestep collapse via correct $\rho \cdot Y$ initialization.
* **CFL Timestep Diagnostics:** Instrumentated `compute_courant_step()` in [HydroSolver.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/HydroSolver.cpp) to report gravity-based vs. velocity-based CFL constraints.

---

## 🚩 Phase 36: Bug Fixes & Release Build Performance Optimization (Completed) 🚀✨
* **Passive Scalar Initializer:** Enabled setting initial passive scalar fractions from config files in [Initializer.cpp](file:///home/bgkang/Projects/RAMSES-CPP/src/Initializer.cpp).
* **Release Build Speedups:** Set default CMAKE_BUILD_TYPE to Release, introducing `-O3` optimizations which yielded a **9.2x speedup** (simulation time for `advect1d` reduced from 46.5s to 5.05s).

---

## 🚩 Phase 35: Dynamic Grid Resizing & Sub-cycling (Completed) 🚀
* **AmrGrid Dynamic Storage:** Implemented [AmrGrid::resize_grids](file:///home/bgkang/Projects/RAMSES-CPP/src/AmrGrid.cpp) to dynamically reallocate and shift per-cell/per-grid arrays, removing hardcoded grid limits and resolving grid allocation overflows during refinement.
* **Dynamic dt Calculation:** Consolidated CFL calculation to determine `nsub1_max_level` dynamically across different level configurations.

---

## 🚩 Phase 34: 0-based Parity & Snapshot Alignment (Completed) 🛰️
* **Tracer Particles:** Added classical tracer particle dynamics with in-place density-proportional initialization.
* **Binary Format Alignment:** Aligned record sizes in `RamsesWriter` for AMR grids, boxlen headers, and hydro levels to matches legacy visualizer expectation.
* **MUSCL Slopes:** Implemented velocity-weighted Superbee/Ultrabee limiters for 1D.

---

## 🚩 Phases 30-33: Unified Indexing & Tree Growth Infrastructure (Completed) 🧬
* **Unified Level Indexing (Phase 33):** Re-indexed coarse level base to 0, ensuring Level 0 represents the coarse grid and Level 1+ represents refined octs.
* **Godunov & Metadata Alignment (Phase 33.1):** Resolved uninitialized memory access at level interfaces. Fixed `header_*.txt` structure to support legacy parsing scripts.
* **Cell Center Coordination (Phase 32):** Fixed coordinate calculations in `get_cell_center` for 1D/2D and established neighbor grid linking in `make_grid_fine`.
* **Testing Infrastructure (Phase 31):** Mapped suite flags to CMake and fixed vector out-of-bounds in `TreeUpdater::smooth_fine`.
* **Modern Infrastructure (Phase 30):** Integrated [Hdf5Writer](file:///home/bgkang/Projects/RAMSES-CPP/src/Hdf5Writer.cpp) and comoving light-cone shell routines.

---

## 🚩 Phases 24-29: Advanced Physics & Star/Sink/Feedback Solvers (Completed) 🌟
* **Clump Finder & Sinks (Phase 28-29):** Ported peak-finding structure analysis and sink particle accretion, angular momentum transfers, and sink merging.
* **Star Formation & Feedback (Phase 26-27):** Implemented [StarSolver](file:///home/bgkang/Projects/RAMSES-CPP/src/StarSolver.cpp) with Poisson-based spawning, gas depletion conservation, supernova feedback energy injection, and metal passive scalar yields.
* **Stability & Physics (Phase 24-25):** Added barotropic EOS models (Isothermal, Polytrope), co-moving self-gravity Poisson calculations, comoving light-cone crossing, and stochastic turbulence mode-sum spectral driving.
* **Relativistic Hydrodynamics (RHD) (Phase 24):** Added relativistic primitive recovery, specialized HLLC/HLL solvers, and 'TM' EOS support.

---

## 🚩 Phases 1-23: Foundation & Legacy Code Migration (Completed) 📂
* **MPI & Load Balancing (Phase 23):** Built Hilbert curve spatial partitioning and asynchronous MPI exchanges for boundaries and particle migration.
* **Gravity (Phase 21):** Implemented Cloud-in-Cell mass assignment, multi-grid Poisson solver, and binary Grafic IC support.
* **Radiative Transfer (Phase 20):** Integrated RT M1 closure transport and chemical non-equilibrium ion fraction tracking.
* **MHD HLLD (Phase 14-19):** Ported HLLD solver, constrained transport magnetic field maintenance, and ISM cooling.
* **Core Base (Phase 1-13):** Designed C++ linked-list octree [AmrGrid](file:///home/bgkang/Projects/RAMSES-CPP/include/ramses/AmrGrid.hpp), MUSCL-Hancock base solver, and legacy snapshot reader.

---
📜 *The journey from 1980s Fortran to 21st-century C++.* 🚀
