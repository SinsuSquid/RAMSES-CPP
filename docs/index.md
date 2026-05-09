---
layout: default
title: Home
---

<img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan.png?raw=true" style="width: 300px; height: auto; float: right; margin-left: 20px;">

# RAMSES-CPP Documentation

Welcome to the documentation for **RAMSES-CPP**, a high-performance, modern C++17 port of the renowned **RAMSES-2025** Adaptive Mesh Refinement (AMR) code.

## Overview

RAMSES-CPP provides a physically consistent, modern C++ alternative to the original Fortran implementation of RAMSES. It maintains **strict bit-perfect parity** in data structures and I/O while offering improved maintainability, modularity, and scalability.

## 🚀 Current Status: Stable Distributed Engine

The project has achieved full maturity in its distributed physics engine. It supports MPI-scaled simulations with full sub-cycling and achievement of bit-perfect snapshots across multiple HPC benchmarks.

### Key Features
- **Multi-Dimensional Build System:** Simultaneously generate `ramses_1d`, `2d`, and `3d` executables from a single CMake configuration. 🛠️
- **Restored AMR Engine:** Constant-time $O(1)$ neighbor lookups and robust tree management with unified 1-based indexing. 🧬
- **Advanced Physics:** Fully ported modules for MUSCL-Hancock Hydro, HLLD/Flux-CT Magnetohydrodynamics (MHD), and M1 Radiative Transfer (RT). 🛰️
- **Distributed Gravity:** Multi-grid Poisson solver coupled with a Cloud-In-Cell (CIC) particle-mesh system for cosmological dynamics. 🌌
- **AI-Assisted Patch System:** Translate legacy Fortran overrides into modern C++ using the `ramses-patch-porter`. 🧠
- **BIT-PERFECT I/O:** Exact binary record compatibility with legacy unformatted Fortran snapshots and existing visualization tools. 📊

## Contents

- [Installation](installation.md) - Build instructions and prerequisites.
- [Usage](usage.md) - How to run simulations and analyze results.
- [Architecture](architecture.md) - Deep dive into the central grid manager and polymorphic solvers.
- [Magnetohydrodynamics (MHD)](mhd.md) - Flux-CT and staggered grid details.
- [Poisson Solver (Gravity)](poisson.md) - Multi-grid comoving self-gravity.

---
Developed with 💖 by Gemini-chan for Senpai. 🚀✨
