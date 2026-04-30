# Poisson Solver (Self-Gravity)

RAMSES-CPP includes a high-performance Multigrid Poisson solver for computing the self-gravitational potential of the matter distribution.

## Algorithm

The solver uses a Geometric Multigrid approach to solve the Poisson equation:
$$\nabla^2 \phi = 4 \pi G \rho$$

### Components:
- **Smoothing:** Iterative Gauss-Seidel smoothing with Red-Black ordering for efficient convergence and easy parallelization.
- **V-Cycle:** Recursive V-cycle strategy.
- **Restriction:** Coarsening the residual from fine to coarse levels.
- **Prolongation:** Interpolating the correction from coarse to fine levels.

## AMR Integration

The Poisson solver is fully integrated with the AMR engine. It solves the equation on the hierarchy of levels, using the potential from coarser levels as boundary conditions for finer levels.

## Configuration

The Poisson solver is included in the default build. To enable gravity in a simulation, set `gravity=.true.` in the namelist.

### Namelist Parameters
- `gravity`: Boolean to enable/disable self-gravity.
- `niter_smooth`: Number of Gauss-Seidel smoothing iterations per level.
- `epsilon_poisson`: Convergence threshold for the iterative solver.
