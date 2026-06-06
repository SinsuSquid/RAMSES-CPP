---
layout: default
title: Home
---

# RAMSES-CPP Documentation 🚀✨

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/SinsuSquid/RAMSES-CPP)

<p align="center">
  <img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan.png?raw=true" width="300" alt="Sia-chan">
</p>

Welcome to the documentation for **RAMSES-CPP**, a high-performance, modern C++17 port of the renowned **RAMSES-2025** Adaptive Mesh Refinement (AMR) astrophysical code.

---

## 🗺️ Documentation Directory

Use the cards below to navigate the core developer and user guides:

| Guide | Description | Target Audience |
| :--- | :--- | :--- |
| **🚀 [Getting Started](installation.md)** | Prerequisites, compiler versions, and build instructions. | Users / Sysadmins |
| **📖 [Simulation Usage](usage.md)** | Namelists configuration, running jobs, and visualization. | Researchers |
| **🌲 [AMR Grid & Memory](amr_grid.md)** | The C++ octree geometry, contiguous vectors, and dynamic resizing. | Core Developers |
| **🔁 [Subcycling & Time Loop](time_integration.md)** | Recursive timestepping, boundaries, restriction, and CFL controls. | Core Developers |
| **🧲 [Physics Solvers](physics_solvers.md)** | Hydro, MHD staggered Flux-CT, Poisson self-gravity, and RT. | Physics developers |
| **🧪 [Testing & Validation](testing_validation.md)** | Differential testing framework and python validation pipeline. | Developers |

---

## 🏆 Modernization Log (Milestones)

The project has achieved **100% architectural and physics parity** with legacy RAMSES-2025. Here is a summary of the porting milestones:

* **Core Foundation:** Established base AMR grid linking, comoving expansion calculations, and Riemann solver class hierarchies.
* **Build & Memory Safety:** Integrated CMake multi-dimensional targets and implemented the dynamic, layout-preserving `resize_grids` algorithm to eliminate grid-overflow crashes.
* **Magnetohydrodynamics (MHD):** Ported staggered face-centered magnetic fields, the Constrained Transport (Flux-CT) divergence cleaning method, and the HLLD Riemann solver to guarantee machine-precision $\nabla \cdot B = 0$.
* **AMR Refinement & Stencil Rules:** Implemented nested grid verification (`ensure_ref_rules`), gradient-based flagging, and non-thermal energy configurations.
* **Limiters & Jump Conditions:** Integrated Toro's exact jump conditions in the HLLC solver and added new slope limiters.
* **Timing Parity & Regression Verification:** Restructured the core integration loop to run via recursive level-by-level `amr_step` subcycling and validated parity against legacy snapshots across all 36 test cases.
