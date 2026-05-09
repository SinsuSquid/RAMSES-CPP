# RAMSES-CPP Architecture

RAMSES-CPP is a modern C++17 port of the RAMSES AMR code, designed with modularity and extensibility as primary goals. The architecture is built around a centralized AMR grid manager and a collection of polymorphic solvers.

## 🏗️ Core Components

### 1. AmrGrid
The `AmrGrid` class manages the distributed octree structure. It handles:
- Memory allocation for cell and grid (oct) data.
- 1-based indexing parity with legacy Fortran for algorithmic consistency.
- Neighbor finding and tree traversal.
- MPI rank mapping and ghost cell synchronization.

### 2. Polymorphic Solvers
All physics modules (`HydroSolver`, `MhdSolver`, `RtSolver`, etc.) inherit from base abstract classes. This allows for:
- **Solver Factory:** Dynamic selection of physics engines (standard vs. experimental) via the `SolverFactory`.
- **Easy Extension:** New physics can be added by implementing a single derived class without touching the core simulation loop.

### 3. Simulation Core
The `Simulation` class orchestrates the time-stepping loop, handling:
- Sub-cycling across AMR levels.
- Global timestep (CFL) reductions.
- Snapshot dumping via the `RamsesWriter`.
- Distributed load balancing using Hilbert curves.

## 🧠 AI-Assisted Patch System
To support the vast ecosystem of legacy RAMSES patches, the project includes the `ramses-patch-porter`. This tool uses LLMs to translate legacy Fortran `.f90` overrides into C++ polymorphic classes that are automatically discovered and compiled into the engine.

## 📊 I/O Layer
The `RamsesWriter` and `RamsesReader` modules are designed for **Bit-Perfect Parity**. They follow the exact binary record sequence used by Fortran unformatted I/O, ensuring full compatibility with existing Python and IDL visualization tools (e.g., `visu_ramses.py`).

## 🛰️ Distributed Scalability
The engine uses a sophisticated `MpiManager` and `LoadBalancer` to distribute octs across nodes. Communication is minimized using asynchronous MPI primitives, and the `TreeUpdater` ensures that the grid hierarchy remains consistent across rank boundaries.

## 🚩 Future Roadmap (Phase 24+)
To complete the full parity with the legacy RAMSES-2025 ecosystem, the following modules are planned for integration:
- **RhdSolver:** Relativistic hydrodynamics implementation with Newton-Raphson primitive recovery.
- **Turbulence Forcing:** Support for external driving terms in the hydro equations via spectral mode-sum driving.
- **Sink Dynamics:** Robust MPI-aware sink particle creation and management using global reductions and broadcasts.
- **ATON Integration:** High-performance GPU radiation transport.

---
🚀 *Engineered for performance and parity.* 🚀
