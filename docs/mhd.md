---
layout: default
title: Magnetohydrodynamics (MHD)
---

# MHD Module

The RAMSES-CPP MHD module provides a robust implementation of the Ideal MHD equations on a distributed AMR grid.

## 🚀 Current Status: Production Ready
The MHD module is fully ported and integrated into the multi-dimensional build system. It supports 1D, 2D, and 3D configurations with machine-precision divergence maintenance. CMake integration (Phase 31) ensures that the `-DRAMSES_USE_MHD` flag correctly activates all magnetic flux components across the automated test suite.

## 🛠️ Key Implementation Details

### 1. Staggered Grid (Flux-CT)
Magnetic field components ($B_x, B_y, B_z$) are stored on cell faces to facilitate the **Constrained Transport (CT)** method. This ensures that $\nabla \cdot B = 0$ is preserved to machine precision throughout the simulation.

### 2. HLLD Riemann Solver
The module uses the HLLD (Harten-Lax-van Leer Discontinuities) solver, which is specifically designed for MHD. It correctly handles the seven waves of the MHD system, providing superior resolution for rotational discontinuities and entropy waves.

### 3. MUSCL Reconstruction
High-order spatial accuracy is achieved using the MUSCL scheme with slope limiters (MinMod, MC, Superbee). The predictor step includes full transverse coupling and source terms.

## 📈 Divergence Monitoring
The code automatically monitors and logs the maximum divergence of the magnetic field at every coarse step.
```text
Step=50 t= 0.005 dt= 0.0001
[MHD] Max DivB: 1.2e-16 (Parity: OK)
```

## 🧪 Verified Benchmarks
- **1D MHD Shock Tube:** Verified wave propagation and jump conditions.
- **2D Orszag-Tang Vortex:** Verified complex shock interactions and $\nabla \cdot B$ maintenance.

---
🚀 *Stability and Parity across all magnetic fields.* 🚀
