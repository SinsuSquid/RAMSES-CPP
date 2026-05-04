# Magnetohydrodynamics (MHD) in RAMSES-CPP

> **Status:** The MHD module is fully integrated with the `AmrGrid` system and is the primary focus of the Phase 2 multi-dimensional expansion. Core physics, Constrained Transport, and adaptive refinement for magnetic gradients are fully operational and verified against standard benchmarks.

The MHD module in RAMSES-CPP provides a robust, divergence-free implementation of ideal magnetohydrodynamics on adaptive grids.

## Key Features

- **Riemann Solvers:** Supports HLLD (standard for MHD) and LLF (Local Lax-Friedrichs) for high robustness.
- **Constrained Transport (CT):** Uses a staggered grid arrangement for magnetic fields to maintain $\nabla \cdot B = 0$ to machine precision.
- **2nd-Order Predictor-Corrector:** Implements a MUSCL-Hancock scheme with a $dt/2$ time prediction step for 2nd-order accuracy in both space and time.
- **Adaptive Magnetic Refinement:** Supports automatic mesh refinement based on magnetic energy gradients (`err_grad_b2`), essential for capturing MHD shocks and vortices.
- **Dimensional Rotation:** Features a robust state-rotation framework in `MhdSolver::cmpflxm` to ensure physical consistency in multi-dimensional sweeps.

## Implementation Details

### Staggered Grid and Variable Layout
Magnetic fields are stored at the faces of the cells (indices 6-8 for left faces, and $nvar-2$ to $nvar$ for right faces), while hydro variables are stored at cell centers.
- $nvar = 11 + nener$ for MHD simulations.
- Fluxes are computed at cell interfaces and correctly back-rotated to the global coordinate system after the Riemann solver step.

### EMF and CT Update
The Electromotive Forces (EMFs) are computed at cell edges using electric field components from the Riemann fluxes. The face-centered magnetic fields are then updated using the curl of these EMFs:
$$ \frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \mathbf{E} $$
This staggered update ensures that the divergence $\nabla \cdot B$ remains zero to machine precision throughout the simulation.

### Verified Benchmarks
- **1D MHD Shock Tube:** Verified correct wave propagation and jump conditions.
- **2D Orszag-Tang Vortex:** Successfully executed high-resolution adaptive runs (~35,000 leaf cells at Level 9). Verified symmetric density evolution and machine-precision divergence maintenance ($\nabla \cdot B \approx 10^{-16}$).

## Configuration

To enable MHD, use the following CMake flag:
```bash
cmake .. -DRAMSES_USE_MHD=ON
```

In the namelist, ensure `hydro=.true.` and appropriate MHD parameters are set in the `&HYDRO_PARAMS` and `&PHYSICS_PARAMS` blocks.
