---
layout: default
title: Installation
---

# Installation

RAMSES-CPP utilizes a modern CMake-based build system to provide a streamlined installation process.

## Prerequisites

Ensure your system meets the following requirements:

- **CMake:** Version 3.15 or higher.
- **C++ Compiler:** C++17 compliant (GCC 9+, Clang 10+, or MSVC 2019+).
- **MPI (Optional):** Required for distributed-memory parallel execution (e.g., OpenMPI or MPICH).
- **Python 3:** Required for automated testing and visualization scripts.

## Building the Project

The build system is designed to generate dimension-specific executables simultaneously.

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/SinsuSquid/RAMSES-CPP.git
   cd RAMSES-CPP
   ```

2. **Configure and Build:**
   ```bash
   mkdir build && cd build
   cmake .. -DMPI=ON -DSOLVER=mhd -DRT=ON
   make -j$(nproc)
   ```

### Core Build Options

| Flag | Description | Default |
|------|-------------|---------|
| `MPI` | Enable MPI-scaled distributed execution. | `OFF` |
| `SOLVER` | Solver type: `hydro`, `mhd`, or `rhd`. | `hydro` |
| `RT` | Enable the Radiation Transport module. | `OFF` |
| `NPRE` | Floating point precision (4 for float, 8 for double). | `8` |
| `NENER` | Number of non-thermal energy variables. | `0` |
| `NPSCAL` | Number of passive scalar variables. | `0` |
| `NMETALS` | Number of metal species (passive scalars). | `0` |
| `NGROUPS` | Number of RT groups. | `0` |
| `NIONS` | Number of RT ions. | `0` |
| `USE_TURB` | Enable FFTW-based turbulence driving. | `OFF` |
| `ATON` | Enable ATON GPU solver. | `OFF` |
| `GRACKLE` | Enable Grackle cooling library. | `OFF` |
| `NVECTOR` | Size of vector cache. | `32` |

## Generated Executables

After a successful build, the `build/` directory will contain:
- `ramses_1d`: Optimized for 1D physics.
- `ramses_2d`: Optimized for 2D physics.
- `ramses_3d`: Full 3D AMR solver.

---
Developed with 💖 by Gemini-chan. 🚀✨
