# RAMSES-2025 C++ Porting History

This file tracks the architectural decisions and progress of the port from Fortran to C++.

## Phase 0: Infrastructure & Foundation

### [2026-04-22] - Project Initialization
- **C++ Environment:** Created `cpp/` directory structure and `CMakeLists.txt`.
- **Standards:** Targeting C++17.
- **Type Mapping:** Created `Types.hpp` for bit-perfect type alignment.
- **Constants:** Ported foundation into `Constants.hpp` and `Parameters.hpp`.

## Phase 4: AMR Tree Logic & Refinement (Current)

### [2026-04-22] - Tree Traversal & Coarse Refinement
- **Neighbor Search:** Implemented `get_nbor_grids` and `get_nbor_cells` in `AmrGrid` using dimension-agnostic lookup tables.
- **Refinement Logic:** Implemented `TreeUpdater::refine_coarse` with periodic boundary support.
- **Verification:** Successfully verified coarse-level oct creation and neighbor connectivity.

## Phase 5: Validation Infrastructure (Current)

### [2026-04-22] - Fortran Binary Bridge
- **RamsesReader:** Implemented a C++ utility to read RAMSES unformatted Fortran binary files.
- **Snapshot Loading:** Developed `load_amr` and `load_hydro` to reconstruct the full state (Grid + Conservative Physics) from Fortran output.
- **Physics Conversion:** Implemented automatic conversion from primitive variables (saved in RAMSES) to conservative variables (used in C++ `uold`).
- **Parity:** This enables bit-for-bit comparison between the legacy Fortran code and the new C++ port.

## Phase 6: Hydro Solver (Current)

### [2026-04-22] - Godunov Solver & 3D Stencil
- **Riemann & MUSCL:** Core physics (LLF, MUSCL-Hancock, Slopes) implemented and verified.
- **Stencil Assembly:** Implemented `gather_stencil` to construct 6x6x6 local cubes for each oct, providing a unified physics buffer.
- **Interpolation:** Added `interpol_hydro` skeleton for coarse-to-fine prolongation at AMR boundaries.
- **Unsplit Logic:** Implemented a functional 1D unsplit skeleton in `godfine1`, demonstrating the full flow from stencil assembly to `unew` update.

### Architectural Decisions
1. **Oct-Centric Processing:** The solver processes data oct-by-oct to maintain cache locality and match RAMSES' vector-sweep philosophy.
2. **Buffer-Driven Physics:** By assembling a 6x6x6 `LocalStencil`, the physics solvers (MUSCL, Riemann) are shielded from the complexities of the AMR tree navigation.

## Final Summary of C++ Port Initialization
This task has successfully initialized the RAMSES-2025 C++ port with the following verified components:
- [x] Modern C++ Build System (CMake)
- [x] 1-Based Indexing Field Wrappers (Fortran compatibility)
- [x] AMR Grid & Linked-List Tree Structures
- [x] Fortran Namelist Parser
- [x] 1D/2D/3D Hilbert Curve Indexing
- [x] Fortran Binary Bridge (RamsesReader)
- [x] Core Hydro Physics (LLF, MUSCL, TVD Slopes)
