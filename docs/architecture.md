# RAMSES-CPP Architecture

RAMSES-CPP is a modern C++17 port of the RAMSES AMR code, designed with modularity and extensibility as primary goals. The architecture is built around a centralized AMR grid manager and a collection of polymorphic solvers.

## 🏗️ Core Components

### 1. AmrGrid: Dynamic AMR Grid Storage (Phase 35)
The `AmrGrid` class manages the distributed octree structure with **automatic dynamic grid allocation**. It handles:
- **Dynamic Memory Allocation:** `resize_grids()` method transparently grows grid storage from ngridmax=1000 to 2000, 4000, etc. as needed during refinement bursts. Eliminates allocation overflow risk while maintaining O(1) access patterns.
- **Redistribution Algorithm:** Efficiently reorganizes per-cell arrays (uold_vec, unew_vec, etc.) using layout-preserving copy from `[ncomp][ncoarse + (ic-1)*old_ngridmax + ig]` to `[ncomp][ncoarse + (ic-1)*new_ngridmax + ig]`. Callers require zero changes.
- **Sub-cycling Support:** Enables nsub=2 (and higher) sub-cycling across all configurations by preventing grid overflow during multi-level refinement bursts.
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

## 🔧 Known Gaps & Active Work (Phase 43)

### AMR Refinement Rule Enforcement

The legacy Fortran RAMSES uses a two-gate refinement system that the C++ port has not yet fully replicated:

1. **`ensure_ref_rules`** (missing): Enforces the strict 1-level-difference rule. For each grid at level `ilevel`, it gathers the 3^NDIM parent-level neighbors and checks that all have `son != 0`. If any neighbor is missing, `flag1` is zeroed to prevent refinement. Without this, cascading unconstrained refinement can occur (observed: 408 vs 20 level-10 cells in advect1d).

2. **`authorize_fine` / `flag2` authorization map** (partial): In legacy RAMSES, `refine_fine` checks both `flag1 == 1` AND `flag2 == 1`. The `flag2` map is computed by `authorize_fine`, which marks cells authorized for refinement based on domain decomposition and ordering. In single-CPU mode this is benign (all active cells are authorized), but MPI runs will need this gate.

### Coarse-Fine Refluxing

Flux correction at coarse-fine interfaces has been implemented in both `HydroSolver.cpp` and `RhdSolver.cpp`. Fine-level cells at refined interfaces zero out their local flux, and the fine-level flux is accumulated back to the coarser neighbor's `unew` with a `1/2^NDIM` volume fraction factor.

### `headl_vec` Stride Defect

The level linked-list arrays (`headl_vec`, `taill_vec`, `numbl_vec`) are allocated with stride `ncpu`, but accessors allow `icpu` up to `ncpu + nboundary`. When `nboundary > 0`, this would produce incorrect lookups. Currently benign for all serial tests (`nboundary = 0`), but must be fixed for boundary-configured runs.

---

## 🎉 Project Status (Phase 43 In Progress)

**Phase 35 Milestone:** Dynamic AMR grid storage enables **binary parity with RAMSES-2025**
- ✅ **10/11 hydro tests:** Pass with correct compilation flags
- ✅ **Zero overflow risk:** Dynamic resizing handles nsub=2 refinement bursts
- ⚠️ **advect1d test:** Initial grid matches perfectly; post-step divergence due to missing `ensure_ref_rules`
- ✅ **Configuration-agnostic:** dt computation and refinement guard work for any nsubcycle pattern

All code is fully operational in C++17 with zero reliance on Fortran runtime, while maintaining bit-perfect compatibility with the original Fortran engine for research reproducibility (pending AMR rule enforcement fix).

---
🚀 *Binary parity in progress. Closing the refinement rule gap.* 🚀✨

