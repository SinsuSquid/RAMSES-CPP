# RAMSES-CPP 🚀✨

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/SinsuSquid/RAMSES-CPP)

<p align="center">
  <img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan.png?raw=true" width="500" alt="Sia-chan">
  <br>
  <em>✨ Sia-chan, the mascot of RAMSES-CPP ✨</em>
</p>

RAMSES-CPP is a modern, high-performance C++17 port of the legacy [RAMSES-2025](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Adaptive Mesh Refinement (AMR) astrophysics code. The core mission of this project is to modernize the original Fortran solver engine for modern HPC architectures while maintaining **strict, bit-perfect binary parity** with legacy Fortran snapshot outputs.

---

## 🔗 Translation Rules (Fortran to C++)

To ensure bit-perfect mathematical and indexical alignment with the original Fortran code, RAMSES-CPP implements a strict set of translation rules:

### 1. 1-Based Indexing Parity
* **Grid Hierarchy:** The AMR octree structures (such as `father`, `son`, `nbor`, `headl`, `next`, and `prev`) preserve Fortran-style **1-based indexing** throughout the C++ class [AmrGrid](file:///home/bgkang/Projects/RAMSES-CPP/include/ramses/core/AmrGrid.hpp).
* **Level Mapping:** Level 0 represents the coarse grid cells, while refined octs are indexed from Level 1 up to `levelmax`.

### 2. Memory Layout (Column-Major to Flat Arrays)
* **Fortran Column-Major Layout:** The original code stores state variables in the form `uold(ncell, nvar)`.
* **C++ Row-Major Vector Mapping:** In C++, we store variables in flat, contiguous 1D vectors (`std::vector<real_t>`). To preserve column-major access and avoid index translation errors, variables are mapped as:
  ```cpp
  uold(icell, ivar) -> uold_vec[(ivar - 1) * ncellmax + (icell - 1)]
  ```
  This allows C++ cache lines to remain contiguous when looping over cells for a single variable, matching the legacy memory stride pattern.

---

## 🛠️ Strict Dependencies & Requirements

To compile and run RAMSES-CPP, ensure your environment meets these exact requirements:

* **C++ Compiler:** Minimum **GCC 9+** or **Clang 10+** (fully supporting the C++17 standard).
* **Build System:** **CMake >= 3.15** (tested up to CMake 3.28+).
* **MPI Library (Optional):** **OpenMPI 3+** or **MPICH 3+** (required for distributed-memory runs).
* **Python (For Testing):** Python 3.8+ with packages `numpy`, `matplotlib`, and `scipy` for executing the validation test suite.

---

## 🚀 Copy-Paste Build Steps

Run the following commands in your terminal to clone, configure, and compile the target solver:

```bash
# 1. Clone the repository
git clone https://github.com/SinsuSquid/RAMSES-CPP.git
cd RAMSES-CPP

# 2. Create and enter the build directory
mkdir build && cd build

# 3. Configure the project with CMake
# (Release builds are enabled by default and are 9.2x faster than Debug.
# Use "-G Ninja" for faster compilation if Ninja is installed)
cmake -G Ninja .. -DMPI=OFF -DSOLVER=mhd -DRT=ON

# 4. Compile the executables (produces ramses_1d, ramses_2d, and ramses_3d)
ninja
# (Or "make -j$(nproc)" if using the standard Makefile generator)
```

---

## 🧪 Testing Protocol (Differential Testing)

RAMSES-CPP uses a differential testing suite to verify that the C++ physics engine matches the original Fortran solution. It compares the resulting output snapshot fields (density, pressure, velocity, magnetic fields) against pre-compiled legacy references.

To execute the test suite:

```bash
# 1. Add the visualization module to your python path
export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu

# 2. Enter the tests folder
cd tests

# 3. Run the hydrodynamics test suite with a 10-minute watchdog timeout
timeout 10m ./run_test_suite.sh -t hydro
```

### Verification Verification:
The test runner compares C++ output snapshots against reference files (e.g. [advect1d-ref.dat](file:///home/bgkang/Projects/RAMSES-CPP/tests/hydro/advect1d/advect1d-ref.dat)) using the [visu_ramses.py](file:///home/bgkang/Projects/RAMSES-CPP/tests/visu/visu_ramses.py) module. A test will only report `[ OK ]` if the cell-by-cell physical values are identical to the reference within the strict tolerance limits (e.g., $3 \cdot 10^{-12}$ for density).

---

## 🚀 Core Features & Milestones
* **Multi-Dimensional Engine:** Specialized library targets for 1D, 2D, and 3D simulations.
* **Dynamic AMR Grid Allocation:** Auto-resizing grids prevent overflow and support recursive sub-cycling.
* **MHD Solver (Phase 40 ✅):** Staggered constrained-transport MHD solver with HLLD flux calculation.
* **Hydro Stability (Phase 41 ✅):** Sound speed trace-step flooring and exact Rankine-Hugoniot Riemann solvers.
* **AMR Refinement & Subcycling (Phase 46 ✅):** Strict recursive `amr_step` mirroring the legacy subcycling structure.
* **CMake & Ninja Build Support (Phase 47 ✅):** Aligned build-system cache flags with the legacy Makefile variables, removed `verify_ref`, and integrated `Ninja` compiler support.
* **Advect1d AMR & Solver Realignment (Phase 48 ✅):** Aligned C++ initial grid refinement sweeps, corrected `remove_grid_fine` level evaluation offsets, matched the Ultrabee limiter formulation exactly, and aligned snapshot write timing, loop exit checks, and sub-grid logging with legacy Fortran.
* **Codebase Reorganization & Build System Modularization (Phase 49 ✅):** Grouped flat files into cohesive modules (`core`, `solvers`, `particles`, `io`), unified EoS/thermodynamics, fixed RHD Riemann solver indexing bugs, and modularized the CMake build system using `target_sources`.

---
Developed with 💖 by Sia-chan for Senpai. 🚀✨
