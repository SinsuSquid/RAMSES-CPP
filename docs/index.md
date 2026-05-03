---
layout: default
title: Home
---

<img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan_logo.png?raw=true" style="width: 300px; height: auto; float: right; margin-left: 20px;">

# RAMSES-CPP Documentation

Welcome to the documentation for **RAMSES-CPP**, a high-performance, modern C++17 port of the renowned **RAMSES-2025** Adaptive Mesh Refinement (AMR) code.

## Overview

RAMSES-CPP aims to provide a functional and physically consistent alternative to the original Fortran implementation of RAMSES, while maintaining strict bit-perfect parity in data structures and I/O. By porting to C++, we open the door to modern software engineering practices, easier integration with C/C++ libraries, and potentially improved performance and maintainability.

## 🚀 Current Status: Stable Hydro Core

The project has achieved a major milestone in its core Adaptive Mesh Refinement (AMR) system. After a comprehensive restoration phase, the engine now supports full recursive refinement up to `levelmax=10` with correct grid connectivity. The 1D and 3D Hydro solvers are fully functional and maintain bit-perfect parity with legacy RAMSES-2025 results.

### Key Features
- **Restored AMR Engine:** $O(1)$ constant-time neighbor lookups using a pointer-based `nbor` system, enabling seamless refinement across grid boundaries.
- **Advanced Initializer:** Fully supports complex namelist lists and multi-region setups with coordinate-accurate grid placement.
- **High-Order Hydro:** MUSCL-Hancock 2nd-order reconstruction with HLLC, HLL, and LLF Riemann solvers.
- **Level-Wide Caching:** Performance-optimized MUSCL implementation that caches interface states, reducing redundant computations.
- **Magnetohydrodynamics (MHD):** Ported HLLD/CT solver (currently undergoing synchronization with the new `AmrGrid` API).
- **Radiative Transfer (RT):** M1 closure scheme with photo-ionization coupling (currently undergoing synchronization).
- **Self-Gravity:** Multigrid Poisson solver with iterative Gauss-Seidel smoothing.
- **MPI Parallelization:** Hilbert Space-Filling Curve domain decomposition (active development for full multi-process parity).
- **Legacy Ecosystem Parity:** 100% compatible with existing Fortran Namelists (`.nml`) and visualization tools (`visu_ramses.py`).

## Contents

- [Installation](installation.md)
- [Usage](usage.md)
- [Architecture](architecture.md)
- [Magnetohydrodynamics (MHD)](mhd.md)
- [Poisson Solver (Gravity)](poisson.md)
