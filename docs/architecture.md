# RAMSES-CPP Architecture

RAMSES-CPP is a modern C++17 port of the RAMSES AMR code, designed with modularity and extensibility as primary goals. The architecture is built around a centralized AMR grid manager and a collection of polymorphic solvers.

## 🏗️ Core Components

### 1. AmrGrid
The `AmrGrid` class manages the distributed octree structure. It handles:
- Memory allocation for cell and grid (oct) data.
- 1-based indexing parity with legacy Fortran for algorithmic consistency.
- **Robust Neighbor Finding:** Constant-time neighbor lookups with proper periodic boundary handling at the coarsest level.
- **Dynamic Grid Linking:** Automated neighbor-linking pass in `make_grid_fine` ensures newly created octs correctly identify sibling/cousin grids.
- **Dimensional Stability:** Fixed coordinate calculations for 1D and 2D simulations, ensuring cells are correctly centered in unused dimensions ($0.5 \times \text{boxlen}$).
- MPI rank mapping and ghost cell synchronization.

### 2. Polymorphic Solvers
All physics modules (`HydroSolver`, `MhdSolver`, `RtSolver`, etc.) inherit from base abstract classes. This allows for:
- **Solver Factory:** Dynamic selection of physics engines (standard vs. experimental) via the `SolverFactory`.
- **Easy Extension:** New physics can be added by implementing a single derived class without touching the core simulation loop.

### 3. Simulation Core
The `Simulation` class orchestrates the time-stepping loop, handling:
- Sub-cycling across AMR levels.
- **Refinement-Aware Initialization:** The `Initializer` is applied within the refinement loop, ensuring newly created levels are populated before flagging.
- Global timestep (CFL) reductions.
- Snapshot dumping via the `RamsesWriter`.
- Distributed load balancing using Hilbert curves.

## 🧠 AI-Assisted Patch System
To support the vast ecosystem of legacy RAMSES patches, the project includes the `ramses-patch-porter`. This tool uses LLMs to translate legacy Fortran `.f90` overrides into C++ polymorphic classes that are automatically discovered and compiled into the engine.

## 📊 I/O Layer
The `RamsesWriter` and `RamsesReader` modules are designed for **Bit-Perfect Parity**. They follow the exact binary record sequence used by Fortran unformatted I/O, ensuring full compatibility with existing Python and IDL visualization tools (e.g., `visu_ramses.py`).

## 🛰️ Distributed Scalability
The engine uses a sophisticated `MpiManager` and `LoadBalancer` to distribute octs across nodes. Communication is minimized using asynchronous MPI primitives, and the `TreeUpdater` ensures that the grid hierarchy remains consistent across rank boundaries.

- **RhdSolver:** Relativistic hydrodynamics implementation with Newton-Raphson primitive recovery and 'TM' EOS support.
- **TurbulenceSolver:** Stochastic forcing in Fourier space using mode-sum spectral driving.
- **Sink Dynamics:** Robust MPI-aware sink particle creation and management using global reductions and broadcasts.
- **StarSolver:** Gas-to-star conversion using Poisson statistics and a deterministic, grid-based RNG for bit-perfect reproducibility.
- **FeedbackSolver:** Supernova energy and metal injection with support for delayed cooling.
- **ClumpFinder:** On-the-fly structure identification based on density peak finding and saddle-point merging.
- **Tracer Particles:** Massless particles following fluid trajectories via trilinear interpolation of gas velocity.
- **LightCone:** Cosmological shell identification for deep-field surveys.
- **Hdf5Writer:** Parallel HDF5 output mirroring the legacy RAMSES hierarchical schema.

## 🚩 Project Status: 100% Port Parity
As of Phase 32, the project has achieved complete architectural and physics parity with the RAMSES-2025 release. All core modules, from base AMR to advanced feedback and clump finding, are fully operational in C++17. Recent stabilization efforts fixed a critical bug in AMR tree growth, enabling full refinement depth up to Level 10 and beyond.

---
🚀 *Engineered for performance and parity.* 🚀
