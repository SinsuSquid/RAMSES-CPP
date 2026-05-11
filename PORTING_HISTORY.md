# 🕰️ RAMSES-CPP Porting History

This document tracks the major milestones and architectural shifts during the migration from legacy RAMSES (Fortran) to the modern C++17 distributed engine.

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

## 🚩 Phase 24: New Physics & Parity (Completed)
- **Relativistic Hydrodynamics (RHD):** Ported `legacy/rhd/` to create a modern `RhdSolver`. Implemented Newton-Raphson primitive recovery and HLLC/HLL/LLF solvers with 'TM' EOS support.
- **Turbulence Driving:** Ported forcing routines from `legacy/turb/`. Implemented `TurbulenceSolver` with Mode-Sum spectral driving as a robust, dependency-free fallback.
- **Sink Particle MPI Fix:** Implemented `SinkSolver` with robust cross-rank synchronization. Added `MPI_Allreduce` for accretion and `MPI_Bcast` for coordinated creation, resolving critical numerical stalls.
- **BIT-PERFECT Alignment (24.1):** Standardized `RamsesWriter` record ordering (`ilevel -> ibound -> ic -> ivar`) and grid coordinate scaling to match legacy RAMSES binary format exactly.

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
