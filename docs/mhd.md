# Magnetohydrodynamics (MHD) in RAMSES-CPP

The MHD module in RAMSES-CPP provides a robust, divergence-free implementation of ideal magnetohydrodynamics on adaptive grids.

## Key Features

- **Riemann Solvers:** Supports HLLD (standard for MHD) and LLF (Local Lax-Friedrichs) for high robustness.
- **Constrained Transport (CT):** Uses a staggered grid arrangement for magnetic fields to maintain $\nabla \cdot B = 0$ to machine precision.
- **MUSCL Reconstruction:** 2nd-order spatial accuracy using MUSCL-Hancock with various TVD slope limiters.
- **Unsplit Integrator:** A fully unsplit 3D integration scheme for improved accuracy and stability.

## Implementation Details

### Staggered Grid
Magnetic fields are stored at the faces of the cells, while hydro variables (density, momentum, energy) are stored at the cell centers. This arrangement is essential for the Constrained Transport method.

### EMF Computation
The Electromotive Forces (EMFs) are computed at the cell edges using Riemann fluxes at the faces. These EMFs are then used to update the face-centered magnetic fields, ensuring that the induction equation is solved in a way that preserves the divergence-free property.

### Divergence Maintenance
The code includes automated monitoring of $\nabla \cdot B$. In multi-dimensional tests like the Orszag-Tang vortex, the divergence is maintained at the level of machine round-off.

## Recent Multi-Dimensional Stability Overhaul

The MHD solver recently underwent a major stabilization phase to resolve numerical explosions in high-gradient multi-dimensional flows.

- **Flux-Consistent EMFs:** Rewrote the EMF construction logic to ensure electric fields are derived from the correct Riemann flux components, maintaining mathematical consistency with the induction equation.
- **Corrected Staggered Signs:** Precisely aligned the signs of the curl operator in the Constrained Transport update with the staggered grid arrangement.
- **MUSCL-Hancock Refinement:** Enhanced the time-prediction step to improve stability at sharp gradients without introducing excessive diffusion.
- **Verified benchmarks:** Successfully executed the 2D Orszag-Tang benchmark through $t=0.5$ at $128^2$ resolution with bit-perfect divergence maintenance.

## Configuration

To enable MHD, use the following CMake flag:
```bash
cmake .. -DRAMSES_USE_MHD=ON
```

In the namelist, ensure `hydro=.true.` and appropriate MHD parameters are set in the `&HYDRO_PARAMS` and `&PHYSICS_PARAMS` blocks.
