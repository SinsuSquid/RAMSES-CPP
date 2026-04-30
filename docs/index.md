---
layout: default
title: Home
---

<img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan_logo.png?raw=true" style="width: 300px; height: auto; float: right; margin-left: 20px;">

# RAMSES-CPP Documentation

Welcome to the documentation for **RAMSES-CPP**, a high-performance, modern C++17 port of the renowned **RAMSES-2025** Adaptive Mesh Refinement (AMR) code.

## Overview

RAMSES-CPP aims to provide a functional and physically consistent alternative to the original Fortran implementation of RAMSES, while maintaining strict bit-perfect parity in data structures and I/O. By porting to C++, we open the door to modern software engineering practices, easier integration with C/C++ libraries, and potentially improved performance and maintainability.

## 🚀 Current Status: Production-Ready Core

The project has reached a production-ready milestone for its core hydrodynamics and gravity solvers. It is capable of executing full 3D simulations (e.g., Sedov Blast) with recursive sub-cycling and high-fidelity physics.

### Key Features
- **Core AMR Engine:** Bit-perfect port of the linked-list octree structure.
- **3D Unsplit Godunov Solver:** Full directional flux integration with strict conservation.
- **Advanced Hydrodynamics:** HLLC and LLF Riemann solvers with MUSCL-Hancock 2nd-order reconstruction and TVD slope limiters (MinMod, MonCen).
- **Magnetohydrodynamics (MHD):** Ported HLLD and LLF Riemann solvers with Flux-CT for $\nabla \cdot B = 0$ maintenance.
- **Self-Gravity:** Multigrid Poisson solver with iterative Gauss-Seidel smoothing.
- **N-Body Dynamics:** Full Particle-Mesh (CIC) mass assignment and advection logic.
- **MPI Parallelization:** Domain decomposition using Hilbert curve repartitioning.
- **Legacy Ecosystem Parity:** 100% compatible with existing Fortran Namelists (`.nml`), test suites, and visualization tools (`visu_ramses.py`).

## Contents

- [Installation](installation.md)
- [Usage](usage.md)
- [Architecture](architecture.md)
- [Magnetohydrodynamics (MHD)](mhd.md)
- [Poisson Solver (Gravity)](poisson.md)
