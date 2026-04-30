# Architecture of RAMSES-CPP

RAMSES-CPP is designed as a modern C++17 object-oriented port of the legacy RAMSES Fortran code. It prioritizes maintainability, type safety, and modularity while maintaining bit-perfect parity with the original simulation logic.

## Core Components

The engine is built around several key classes that manage the simulation lifecycle, the AMR grid, and the physics solvers.

### 1. `Simulation` (The Orchestrator)
The `Simulation` class is the central driver. It manages the global time-stepping, level sub-cycling, and coordinates between the grid management and the physics solvers.
- **Path:** `src/Simulation.cpp`, `include/ramses/Simulation.hpp`
- **Key Responsibilities:**
    - Parsing namelists (via `Parameters`).
    - Initializing the grid (via `Initializer`).
    - Running the main simulation loop.
    - Managing recursive level-stepping.

### 2. `AmrGrid` (The Data Structure)
`AmrGrid` encapsulates the linked-list octree structure. In RAMSES-CPP, this is implemented as a collection of `Oct` and `Grid` structures, mimicking the Fortran memory layout but with C++ container management where appropriate.
- **Level Indexing:** The grid uses **0-based level indexing** (Level 0 = Coarse Grid), ensuring consistency across all solvers.
- **Path:** `src/AmrGrid.cpp`, `include/ramses/AmrGrid.hpp`
- **Key Responsibilities:**
    - Storing cell and grid data.
    - Neighbor lookups and tree traversal.
    - Handling periodic and physical boundaries.

### 3. `TreeUpdater` (Grid Evolution)
`TreeUpdater` handles the dynamic nature of AMR: refinement and de-refinement.
- **Path:** `src/TreeUpdater.cpp`, `include/ramses/TreeUpdater.hpp`
- **Key Responsibilities:**
    - Refining cells based on gradients or user-defined criteria.
    - Smoothing the refinement map (ensuring no more than 2:1 level jumps).
    - Updating the octree topology.

### 4. `HydroSolver` and `MhdSolver` (The Physics)
These classes implement the Godunov-type solvers for fluid dynamics.
- **Path:** `src/HydroSolver.cpp`, `src/MhdSolver.cpp`
- **Key Responsibilities:**
    - Stencil gathering (6x6x6 local cubes).
    - MUSCL-Hancock reconstruction.
    - Riemann solver execution (HLLC, LLF, HLLD).
    - Flux integration and cell state updates.
    - (MHD) Constrained Transport (CT) for $\nabla \cdot B = 0$.

### 5. `PoissonSolver` (Gravity)
Implements self-gravity using a Multigrid approach.
- **Path:** `src/PoissonSolver.cpp`, `include/ramses/PoissonSolver.hpp`
- **Key Responsibilities:**
    - Iterative Gauss-Seidel smoothing with Red-Black ordering.
    - V-Cycle execution (Restriction/Prolongation).

### 6. `CoolingSolver` (Thermal Physics)
Handles gas cooling and heating processes.
- **Path:** `src/CoolingSolver.cpp`, `include/ramses/CoolingSolver.hpp`
- **Key Responsibilities:**
    - Conversion between code units and physical cgs units.
    - Implementation of analytic and tabular cooling models (e.g., ISM model).
    - Semi-implicit integration of the stiff energy equation.

### 7. `RtSolver` (Radiation)
Implements Radiative Transfer using the M1 closure scheme.
- **Path:** `src/RtSolver.cpp`, `include/ramses/RtSolver.hpp`
- **Key Responsibilities:**
    - Advection of photon number density and fluxes for multiple groups.
    - Interpolation of HLL eigenvalues from pre-computed tables.
    - **Ionization Coupling:** Utilizes `RtChemistry` to solve for Hydrogen/Helium ionization fractions and gas energy feedback.

### 8. `MpiManager` and `LoadBalancer` (Parallelism)
- **`MpiManager`**: A singleton that handles MPI initialization and provides wrappers for common collective operations.
- **`LoadBalancer`**: Implements domain decomposition using the Hilbert Space-Filling Curve. It ensures spatial locality and provides full physical state migration (including Hydro/MHD variables) between MPI ranks.

## Data Flow

1. **Initialization:** `Simulation` reads parameters and calls `Initializer` to set up the initial conditions on the `AmrGrid`.
2. **Main Loop:** 
    - `Simulation` determines the time step `dt`.
    - Physics solvers (`HydroSolver`/`MhdSolver`) update the cell states.
    - `TreeUpdater` checks for refinement criteria and updates the grid structure.
    - `PoissonSolver` (if enabled) computes the gravitational potential.
    - `RamsesWriter` periodically saves the state to disk.
3. **Sub-cycling:** The code recursively steps through AMR levels, ensuring that finer levels are updated more frequently with smaller time steps.

## I/O Parity

A critical feature of RAMSES-CPP is its ability to read and write unformatted Fortran binaries. This is handled by `RamsesReader` and `RamsesWriter`, ensuring that C++ snapshots are identical to Fortran snapshots, allowing for seamless use of the existing RAMSES ecosystem.
