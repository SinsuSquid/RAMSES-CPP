---
layout: default
title: Installation
---

# Installation Guide

RAMSES-CPP utilizes a modern, platform-independent CMake build system.

---

## 🛠️ Prerequisites & Requirements

Before building the code, ensure your environment meets the following requirements:

* **C++ Compiler:** Minimum **GCC 9+** or **Clang 10+** (supporting the C++17 standard).
* **Build System:** **CMake >= 3.15** (tested up to CMake 3.28+).
* **MPI Library (Optional):** Required for parallel distributed runs (e.g., OpenMPI 3+ or MPICH 3+).
* **Python (For Testing):** Python 3.8+ with packages `numpy`, `matplotlib`, and `scipy` for executing validation plots.

---

## 🚀 Step-by-Step Compilation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/SinsuSquid/RAMSES-CPP.git
   cd RAMSES-CPP
   ```

2. **Configure with CMake:**
   Initialize a build folder and configure your compilation flags. By default, RAMSES-CPP performs an optimized **Release build** (which compiles with `-O3` and is **9.2x faster** than a Debug build). You can append `-G Ninja` if you have the `ninja` build tool installed to speed up compiling/re-compiling:
   ```bash
   mkdir build && cd build
   cmake -G Ninja .. -DMPI=OFF -DSOLVER=mhd -DRT=ON
   ```

3. **Compile:**
   Compile the executables. The build system will generate the solver customized for the dimension specified in `NDIM` (defaults to 3D):
   ```bash
   ninja
   # (Or "make -j$(nproc)" if configured without the Ninja generator)
   ```

---

## ⚙️ Compilation Flags Mapped to Legacy RAMSES

To maintain compatibility with legacy build pipelines, the CMake variables are named identically to the original Fortran `Makefile` parameters:

| Flag | Description | Values | Default |
| :--- | :--- | :--- | :--- |
| `NDIM` | Dimensionality of the solver. | `1`, `2`, `3` | `3` |
| `NPRE` | Real precision size. | `4` (float), `8` (double) | `8` |
| `SOLVER` | Physics solver engine type. | `hydro`, `mhd`, `rhd` | `hydro` |
| `RT` | Enable M1 Radiation Transport solver. | `ON`, `OFF` | `OFF` |
| `MPI` | Enable MPI-scaled parallel execution. | `ON`, `OFF` | `OFF` |
| `USE_TURB` | Enable FFTW-based turbulence driving. | `ON`, `OFF` | `OFF` |
| `ATON` | Enable ATON GPU radiative transfer. | `ON`, `OFF` | `OFF` |
| `GRACKLE` | Enable Grackle chemistry/cooling library. | `ON`, `OFF` | `OFF` |
| `NENER` | Number of additional non-thermal energies. | Integer $\ge 0$ | `0` |
| `NPSCAL` | Number of passive scalar variables. | Integer $\ge 0$ | `0` |
| `NMETALS` | Number of metal variables (adds to passive scalars). | Integer $\ge 0$ | `0` |
| `NGROUPS` | Number of RT radiation energy groups. | Integer $\ge 0$ | `0` |
| `NIONS` | Number of RT chemical species ions. | Integer $\ge 0$ | `0` |
| `NVECTOR` | Size of vector cache. | Integer $\ge 0$ | `32` |
| `LONGINT` | Use 64-bit integers for particle/grid IDs. | `ON`, `OFF` | `OFF` |
