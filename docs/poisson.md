---
layout: default
title: Poisson Solver
---

# Self-Gravity Module

RAMSES-CPP includes a specialized module for calculating self-gravity on Adaptive Mesh Refinement (AMR) grids.

## 🚀 Current Status: Distributed Scaling
The Poisson solver is fully functional and supports distributed-memory parallel execution via MPI. It is tightly coupled with the `ParticleSolver` for cosmological simulations.

## 🛠️ Key Implementation Details

### 1. Multi-Grid Method
The solver utilizes a geometric multi-grid algorithm to solve the comoving Poisson equation.
- **Smoother:** Gauss-Seidel iterations with Red-Black ordering for efficient convergence.
- **V-Cycle:** Recursive restriction to coarser levels and prolongation (interpolation) back to finer levels.

### 2. Comoving Source Scaling
For cosmological runs, the solver correctly scales the density source term using the expansion factor $a(t)$:
$$\nabla^2 \phi = 1.5 \Omega_m a \frac{\rho - \bar{\rho}}{\rho}$$

### 3. Particle-Mesh (PM) Coupling
Gravity is coupled to the N-body particle system using the **Cloud-In-Cell (CIC)** projection.
- **Mass Assignment:** Particle masses are projected onto the grid to form the density source.
- **Force Interpolation:** The gravitational potential gradient (force) is interpolated back to particle positions for the next "move" step.

## 📈 Performance
The multi-grid solver typically achieves convergence within 10-15 V-cycles for most benchmarks, even with deep AMR refinement.

---
🚀 *Scaling gravity across the cosmos.* 🚀
