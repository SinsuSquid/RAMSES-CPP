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
   cmake .. -DRAMSES_USE_MPI=ON -DRAMSES_USE_MHD=ON -DRAMSES_USE_RT=ON
   make -j$(nproc)
   ```

### Core Build Options

| Flag | Description | Default |
|------|-------------|---------|
| `RAMSES_USE_MPI` | Enable MPI-scaled distributed execution. | `OFF` (Auto-detected) |
| `RAMSES_USE_MHD` | Enable the Magnetohydrodynamics module. | `OFF` |
| `RAMSES_USE_RT` | Enable the Radiation Transport module. | `OFF` |
| `RAMSES_PRECISION` | Floating point precision (4 for float, 8 for double). | `8` |
| `RAMSES_NENER` | Number of non-thermal energy variables. | `0` |
| `RAMSES_NPSCAL` | Number of passive scalar variables. | `0` |

## Generated Executables

After a successful build, the `build/` directory will contain:
- `ramses_1d`: Optimized for 1D physics.
- `ramses_2d`: Optimized for 2D physics.
- `ramses_3d`: Full 3D AMR solver.
- `verify_ref`: Parity verification tool.

---
Developed with 💖 by Gemini-chan. 🚀✨
