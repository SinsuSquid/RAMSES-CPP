---
layout: default
title: Installation
---

# Installation

Building RAMSES-CPP is designed to be straightforward using modern CMake.

## Prerequisites

To build RAMSES-CPP from source, ensure you have the following installed on your system:

- **CMake** (version 3.15 or higher)
- **C++17 compliant compiler** (e.g., GCC 9+, Clang 10+, Apple Clang)
- **MPI (Message Passing Interface)** (Optional, but recommended for parallel execution. e.g., OpenMPI, MPICH)

## Building from Source

1. **Clone the repository:**
   ```bash
   git clone https://github.com/SinsuSquid/RAMSES_CPP.git
   cd RAMSES_CPP
   ```

2. **Create a build directory and configure with CMake:**
   ```bash
   mkdir build && cd build
   cmake ..
   ```

3. **Compile the project:**
   ```bash
   make -j$(nproc)
   ```

   This will generate the dimensional executables (e.g., `ramses_1d`, `ramses_3d`) and the reference verification tool `verify_ref` in the `build/` directory.

## CMake Configuration Options

You can customize the build using the following CMake flags:

- `-DRAMSES_NDIM=[1|2|3]`: Set the number of dimensions (default is 3).
- `-DRAMSES_NENER=N`: Set the number of non-thermal energy variables (default is 0).
- `-DRAMSES_USE_MHD=[ON|OFF]`: Enable or disable the Magnetohydrodynamics module (default is OFF).
- `-DRAMSES_USE_MPI=[ON|OFF]`: Force enable/disable MPI support (usually auto-detected).

Example: Building for 2D MHD
```bash
cmake .. -DRAMSES_NDIM=2 -DRAMSES_USE_MHD=ON
make -j
```

## Parallel Support (MPI)

If an MPI implementation is installed and detected by CMake, the build system will automatically define `RAMSES_USE_MPI`. This enables the `LoadBalancer` class for distributed execution and compiles the code to run across multiple MPI ranks. If MPI is not found, the code safely falls back to a sequential build.
