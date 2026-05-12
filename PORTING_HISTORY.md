# 🕰️ RAMSES-CPP Porting History

This document tracks the major milestones and architectural shifts during the migration from legacy RAMSES (Fortran) to the modern C++17 distributed engine.

## 🚩 Phase 34: Tracer Particle Implementation (Active) 🛰️
- **Tracer Physics Engine:** Implemented classical tracer particles that follow gas velocity using trilinear interpolation from the AMR grid.
- **In-place Initialization:** Ported the `load_tracers_inplace` logic, allowing tracers to be spawned proportionally to gas density at simulation start.
- **Particle Metadata Parity:** Expanded the `ParticlePacket` and type system to include `family` and `tag` fields, ensuring tracers are correctly identified and preserved across MPI rank boundaries.
- **I/O Fixes:** Corrected `RamsesWriter` to properly handle sparse particle arrays by only writing active particles (based on `idp > 0`), maintaining binary compatibility with legacy tools.
- **C++17 Optimization:** Utilized `<random>` for deterministic, seed-based tracer spawning, replacing the legacy Fortran RNG while maintaining reproducible spatial distributions.

## 🚩 Phase 33.1: Stability & I/O Parity (Completed) 🛠️
- **Godunov Solver Stabilization:** Fixed uninitialized memory access in `godunov_fine` where interface cells between levels were reading garbage traces. Implemented robust fallbacks to raw cell primitives at AMR interfaces, resolving simulation stalls and tiny `dt` issues.
- **Snapshot Metadata Parity:** Ensured `header_*.txt` and file descriptors (`hydro_file_descriptor.txt`, `part_file_descriptor.txt`) are always produced with correct naming and directory placement, satisfying legacy visualization and analysis scripts.
- **Fixed Snapshot Mapping:** Corrected the mapping between C++ Level 0 (coarse) and file Level 1, ensuring bit-perfect alignment with legacy snapshot formats and resolving `struct.unpack` buffer errors in `visu_ramses.py`.
- **Recursion Safety:** Fixed a segfault in `Simulation::amr_step` by adding strict boundary checks to prevent out-of-bounds recursion beyond `nlevelmax`.

## 🚩 Phase 33: Unified Level Indexing & Physics Alignment (Completed) 🧬
- **Coarse Level Re-indexing:** Shifted the base AMR level from 1 to 0 to perfectly match legacy RAMSES conventions. Coarse cells (ncoarse) are now Level 0, and the first refined grids are Level 1.
- **Global Refactor:** Updated `TreeUpdater`, `Simulation`, and all physics solvers (`Hydro`, `Mhd`, `Rhd`, `Turbulence`, `Sink`, `Initializer`) to respect the new 0-based level convention.
- **Correction of $dx$ Scaling:** Fixed the cell size calculation to $dx = \text{boxlen} / (nx \cdot 2^{ilevel})$, ensuring bit-perfect parity with Fortran's geometric scaling.
- **Improved Simulation Loop:** Corrected the sub-cycling and density aggregation loops to include Level 0, ensuring mass conservation and Poisson stability across the entire AMR tree.
- **Robust Tree Growth:** Aligned `flag_fine` and `make_grid_fine` logic with the new indexing, enabling deeper and more stable refinement (Level 0 up to `nlevelmax`).

## 🚩 Phase 32: Cell Center Fix & Full AMR Tree Growth (Completed) 🎯
- **Critical Bug: `get_cell_center` Dimension Fix:** Discovered that `AmrGrid::get_cell_center` was computing garbage y/z coordinates for 1D (and z for 2D) simulations by reading uninitialized `xg` entries for unused dimensions. The fix defaults inactive dimensions to `0.5 * boxlen`, which is the center of the unit box. This was the root cause of broken initial conditions — all cells appeared to be inside the high-density region, preventing gradient-based refinement from triggering.
- **Coarse Neighbor Logic:** Replaced the placeholder `get_nbor_cells_coarse` (which returned boundary markers) with proper periodic neighbor indexing using coordinate math, enabling cross-cell gradient computation at the coarsest level.
- **Neighbor Grid Linking in `make_grid_fine`:** Added a full neighbor-linking pass after grid creation, so newly created octs at level N+1 can correctly identify their sibling and cousin grids for ghost zone exchanges and gradient computation.
- **Initializer Loop Fix:** Moved `initializer_->apply_all()` inside the per-level refinement loop so that each newly created level gets properly initialized before gradient flagging at the next level.
- **`cpu_map` Global Init:** Changed default `cpu_map` initialization from 0 to `myid` so all cells are visible to the output reader regardless of when they were created.
- **Result:** AMR tree now grows from level 2 to level 10 (matching Fortran), initial conditions correctly show the density step function, and gradient-based refinement properly cascades through all levels.

## 🚩 Phase 31: Testing Infrastructure Integration & Segfault Fix (Completed)
- **Test Suite CMake Integration:** Updated `tests/run_test_suite.sh` to natively compile RAMSES-CPP using CMake, mapping Fortran flags (`EXEC`, `NDIM`, `SOLVER`) to CMake flags (`-DRAMSES_NDIM`, `-DRAMSES_USE_MHD`).
- **Rule Relaxation:** Updated `GEMINI.md` to explicitly allow modifying `.sh` files inside `tests/` based on user approval.
- **TreeUpdater Segfault Fix:** Diagnosed and fixed a segmentation fault in `TreeUpdater::smooth_fine` where `nbor` array values (which represent neighbor cell indices) were incorrectly treated as grid indices, causing `std::vector` out-of-bounds assertions. Refactored the neighbor lookup logic to properly use `AmrGrid::get_nbor_cells`.
- **Print Optimization:** Fixed a console spam issue in `HydroSolver::hydro_step` by adding a static print guard for `slope_type`.

## 🚩 Phase 30: Modern Infrastructure (Completed) 📂
- **HDF5 Snapshot Suite:** Implemented `Hdf5Writer` with support for hierarchical data groups matching legacy RAMSES schema.
- **Cosmological Light-cone:** Ported light-cone shell identification and particle crossing logic from `legacy/amr/light_cone.f90`.

## 🚩 Phase 29: Clump Finder (Completed) 🏔️
- **Structure Identification:** Ported the iterative peak-finding and saddle-point merging algorithm from `legacy/pm/clump_finder.f90`.
- **Property Computation:** Implemented mass, COM, and velocity dispersion calculations for identified clumps.

## 🚩 Phase 28: Advanced Sink Dynamics (Completed) 🕳️
- **Coordinated Creation:** Implemented the `clump_finder`-linked sink creation logic.
- **Sink Merging:** Added robust sink-sink merging and angular momentum transfer during accretion.

## 🚩 Phase 27: Feedback & Metals (Completed) 💥
- **Supernova Feedback:** Implemented thermal energy injection and delayed cooling (Phase 25) logic.
- **Metal Enrichment:** Added `imetal` passive scalar field and yield-based metal injection from star particles.

## 🚩 Phase 26: Star Formation (StarSolver) (Completed) 🌟
- **Star Spawning Engine:** Implemented `StarSolver` with Poisson-based particle creation logic from `legacy/pm/star_formation.f90`.
- **Deterministic RNG:** Integrated a grid-based spatial hashing seed for `std::mt19937`, ensuring bit-perfect reproducibility across any MPI rank configuration.
- **Gas Depletion:** Implemented mass, momentum, and energy conservation during star spawning, with a 90% gas depletion safety cap.
- **Simulation Integration:** Integrated `StarSolver` into `amr_step` with namelist control (`run_params: star`).

## 🚩 Phase 25: Stability & EOS Refinement
- **Barotropic EOS Suite:** Implemented full support for Isothermal, Polytrope, and Double Polytrope EOS models. Pressure is now correctly recovered from density when `barotropic_eos` is enabled.
- **Poisson Stability Fix:** Integrated the Poisson solve into level-specific `amr_step` to match legacy timing. Added mean density subtraction for stable closed-box simulations.
- **Linear Potential Prolongation:** Replaced naive prolongation with linear interpolation using parent forces, preventing exponential potential growth.
- **Unified Courant Step:** Consolidated global timestep logic into a barotropic-aware `HydroSolver` method, eliminating redundancy and instability.
- **Infrastructure Overhaul:** Fixed `multi_gcov_aggregator.py` for multi-dimensional builds and integrated `timeout` safety into `run_test_suite.sh`.
- **Dynamic AMR Refinement (Phase 25) (Completed) 🔄:** Ported runtime grid adaptation logic (`flag_fine`, `make_grid_fine`, `remove_grid_fine`) into `Simulation::amr_step`. Implemented 1D gradient flagging and `smooth_fine` expansion buffers, enabling the simulation to track moving features (like the advection pulse) at maximum resolution dynamically. ✨💖🚀

## 🚩 Phase 24: New Physics & Parity (Completed)
- **Relativistic Hydrodynamics (RHD):** Ported `legacy/rhd/` to create a modern `RhdSolver`. Implemented Newton-Raphson primitive recovery and HLLC/HLL/LLF solvers with 'TM' EOS support.
- **Turbulence Driving:** Ported forcing routines from `legacy/turb/`. Implemented `TurbulenceSolver` with Mode-Sum spectral driving as a robust, dependency-free fallback.
- **Sink Particle MPI Fix:** Implemented `SinkSolver` with robust cross-rank synchronization. Added `MPI_Allreduce` for accretion and `MPI_Bcast` for coordinated creation, resolving critical numerical stalls.
- **BIT-PERFECT Alignment (24.1) (Completed):** Standardized `RamsesWriter` record ordering (`ilevel -> ibound -> ic -> ivar`) and grid coordinate scaling to match legacy RAMSES binary format exactly. Implemented Superbee slope limiter and fixed coordinate offsets in `visu_ramses` interaction. **Standardized STDOUT** to include legacy-style mesh reports, step timing, and verbose level tracking.

## 🚩 Phase 23: Final Optimization and Parity
- **MPI Manager:** Implemented centralized rank management and asynchronous buffer swaps.
- **Global Reductions:** Integrated MPI-aware CFL timestep calculations and total mass/density reductions.
- **Hilbert Load Balancing:** Ported the Hilbert curve partitioning logic to handle massive octree distributions.
- **Particle Migration:** Implemented asynchronous MPI exchange for particles crossing rank boundaries.

## 🚩 Phase 21: Gravity & Particle Dynamics
- **CIC Projections:** Implemented Cloud-In-Cell mass assignment and force interpolation.
- **Poisson Solver:** Ported the iterative multi-grid Poisson solver for comoving gravitational potential.
- **Grafic ICs:** Added unformatted binary support for Grafic initial conditions (velocities and displacements).
- **Friedman Solver:** Implemented expansion factor tables and growth factor calculations in `Cosmology.cpp`.

## 🚩 Phase 20: Radiation Transport (RT) & Chemistry
- **M1 Closure:** Implemented the M1 moment closure for anisotropic radiation fields.
- **HLL Riemann Solver:** Added a specialized Riemann solver for photon flux.
- **Non-Equilibrium Chemistry:** Ported the `RtChemistry` module for ion fraction tracking (HII, HeII, HeIII).

## 🚩 Phase 14-19: MHD & Stability
- **HLLD Riemann Solver:** Integrated magnetic-aware flux calculations.
- **Constrained Transport (CT):** Implemented $\nabla \cdot B = 0$ maintenance on staggered grids.
- **Stencil Robustness:** Optimized 6x6x6 stencil gathering for high-gradient shocks.
- **ISM Cooling:** Ported Hennebelle (2005) analytic cooling/heating models.

## 🚩 Phase 1-13: Foundations
- **AMR Tree Core:** Initial C++ implementation of the linked-list based octree.
- **Base Hydro:** MUSCL-Hancock scheme with HLLC/LLF Riemann solvers.
- **RamsesReader:** C++ bridge for loading legacy Fortran snapshots.
- **Test Infrastructure:** Initial integration with `visu_ramses.py` and automated suites.

---
📜 *Detailed records recovered and merged by Gemini-chan.* 💖
🚀 *The journey from 1980s Fortran to 21st-century C++.* 🚀
