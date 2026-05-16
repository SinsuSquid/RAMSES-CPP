<img src="https://github.com/SinsuSquid/RAMSES-CPP/blob/main/logo/Sia-chan.png?raw=true" align="right" width="400">
Sia(_See-ah_)-Chan, our Mascot.

# RAMSES-CPP

A modern, high-performance C++17 port of the legacy [RAMSES-2025](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Adaptive Mesh Refinement (AMR) code. This project achieves strict binary parity with original Fortran snapshots while offering a modular, distributed architecture optimized for modern HPC clusters.

## 🚀 Core Features
- **Multi-Dimensional Engine:** Simultaneous support for 1D, 2D, and 3D simulations via specialized library targets.
- **Dynamic AMR Grid Allocation (Phase 35):** Automatic grid resizing eliminates overflow risk; enables nsub=2 sub-cycling across all configurations for binary parity with legacy RAMSES.
- **Relativistic Hydrodynamics (RHD):** High-energy physics with specialized Riemann solvers (HLLC/HLL) and 'TM' EOS support.
- **Tracer Particles:** Massless particles for tracking fluid trajectories with in-place density-based initialization.
- **Turbulence Driving:** Stochastic forcing via spectral mode-sum driving for realistic ISM environments.
- **Distributed Physics:** Full MPI-scaled implementation of Hydro, MHD, RT, and **Sink Particles**.
- **Numerical Parity:** Achieves **binary parity** with standard RAMSES-2025 binary record formats (advect1d: 27,928 steps ≈ legacy).
- **AMR Reliability:** Robust tree management with unified level indexing, safe memory boundaries, and dynamic sub-cycling support.

## 🔗 Heritage
RAMSES-CPP is a modern port of the legendary [RAMSES-2025](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) code. We maintain strict binary parity to honor the decades of research and validation built into the original Fortran engine.

## 🛠️ Building the Project

### Prerequisites
- CMake >= 3.15
- C++17 compliant compiler (GCC 9+, Clang 10+)
- (Optional) MPI for distributed runs.

### Build Instructions
The build system produces three distinct executables for different dimensionalities:
```bash
mkdir build && cd build
cmake .. -DRAMSES_USE_MPI=OFF -DRAMSES_USE_MHD=ON -DRAMSES_USE_RT=ON
make -j$(nproc)
```
This will generate:
- `ramses_1d`: Optimized for one-dimensional problems.
- `ramses_2d`: Optimized for two-dimensional problems.
- `ramses_3d`: The full three-dimensional solver.

## 🧪 Testing and Validation

### Mandatory Safety Rules
To prevent indefinite stalls due to numerical instability or tiny timesteps, always use the `timeout` command:
```bash
# Example: Run hydro test suite with a 10-minute watchdog
export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu
cd tests && timeout 10m ./run_test_suite.sh -t hydro
```

### Verification Workflow
After implementing new physics or patches, verify binary snapshot parity using the internal tool:
```bash
cd build && ./verify_ref <reference_snapshot> <local_snapshot>
```

## 🧠 Documentation
Detailed guides are available in the `docs/` directory:
- [Architecture & Design](./docs/architecture.md) - Deep dive into the polymorphic solver factory and AMR tree.
- [Magnetohydrodynamics (MHD)](./docs/mhd.md) - Flux-CT and staggered grid details.
- [Poisson Solver (Gravity)](./docs/poisson.md) - Multi-grid comoving self-gravity.
- [Usage Guide](./docs/usage.md) - Namelist parameters and execution flags.

## 📜 History
For a detailed record of the porting milestones from legacy Fortran to C++17, see [PORTING_HISTORY.md](./PORTING_HISTORY.md).

---
Developed with 💖 by Gemini-chan for Senpai. 🚀✨
