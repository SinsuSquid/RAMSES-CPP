# 🕰️ RAMSES-CPP Porting History

This document tracks the major milestones and architectural shifts during the migration from legacy RAMSES (Fortran) to the modern C++17 distributed engine.

## 🚩 Phase 23: Final Optimization and Parity (Current)
- **Multi-Dimensional Libraries:** Refactored CMake to build `ramses_lib_1d`, `2d`, and `3d` simultaneously, ensuring macro consistency across the entire call stack.
- **BIT-PERFECT I/O:** Fully aligned `RamsesWriter` binary records with Fortran unformatted I/O. Corrected record padding and metadata offsets for full compatibility with `visu_ramses.py`.
- **Numerical Stability:** Fixed a critical "Ghost Grid" bug where restricted child data was being overwritten by stale parent buffers.
- **Watchdog Mandate:** Integrated `timeout` safety rules into the testing workflow to handle numerical stalls.

## 🚩 Phase 22: Distributed Cosmology & MPI Scaling
- **MPI Manager:** Implemented centralized rank management and asynchronous buffer swaps.
- **Global Reductions:** Integrated MPI-aware CFL timestep calculations and total mass/density reductions.
- **Hilbert Load Balancing:** Ported the Hilbert curve partitioning logic to handle massive octree distributions.

## 🚩 Phase 21: Gravity & Particle Dynamics
- **CIC Projections:** Implemented Cloud-In-Cell mass assignment and force interpolation.
- **Poisson Solver:** Ported the iterative multi-grid Poisson solver for comoving gravitational potential.
- **Grafic ICs:** Added unformatted binary support for Grafic initial conditions (velocities and displacements).

## 🚩 Phase 20: Radiation Transport (RT) & Chemistry
- **M1 Closure:** Implemented the M1 moment closure for anisotropic radiation fields.
- **HLL Riemann Solver:** Added a specialized Riemann solver for photon flux.
- **Non-Equilibrium Chemistry:** Ported the `RtChemistry` module for ion fraction tracking (HII, HeII, HeIII).

## 🚩 Phase 1-19: Foundations
- **AMR Tree Core:** Initial C++ implementation of the linked-list based octree.
- **Base Hydro:** MUSCL-Hancock scheme with HLLC/LLF Riemann solvers.
- **MHD:** Dedner cleaning and CT implementation for magnetic fields.
- **Solver Factory:** Creation of the abstract factory pattern for modular physics injection.

---
🚀 *The journey from 1980s Fortran to 21st-century C++.* 🚀
