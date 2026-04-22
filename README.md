# RAMSES-CPP

A high-performance, modern C++17 port of the **RAMSES-2025** Adaptive Mesh Refinement (AMR) code. 

RAMSES-CPP provides a functional, physically consistent alternative to the original Fortran implementation while maintaining strict binary compatibility for snapshots and restart files.

## 🚀 Status: Production-Ready
The code now supports:
- **3D Unsplit Godunov Hydrodynamics:** HLLC/LLF Riemann solvers with MUSCL-Hancock reconstruction.
- **Multigrid Poisson Solver:** Iterative Gauss-Seidel for self-gravity.
- **N-Body Dynamics:** Full Particle-Mesh (CIC) assignment and advection.
- **AMR Engine:** Recursive sub-cycling with dynamic 3D tree traversal.
- **Ecosystem Parity:** Compatible with legacy RAMSES tools (`visu_ramses.py`) and standard `.nml` namelists.

---

## 🛠 Building the Project

### Prerequisites:
- CMake (>= 3.15)
- C++17 compliant compiler (GCC 9+, Clang, etc.)
- MPI (Optional: detect automatically, defaults to sequential if absent)

### Build Instructions:

1. **Clean build from source:**
```bash
# Create and enter build directory
mkdir build && cd build

# Configure and compile
cmake ..
make -j$(nproc)
```

2. **Parallel Build:**
If you have MPI installed, CMake will automatically detect it and define `RAMSES_USE_MPI`, enabling the `LoadBalancer` for parallel execution.

---

## 🏃 Running a Simulation

Use the standard RAMSES namelist format. The executable is located in your build directory:

```bash
# Example Sedov 3D test
./ramses_main ../namelist/sedov3d.nml
```

Results are saved to `output_XXXXX/`, which is fully compatible with legacy RAMSES visualization tools.

---

## 🗺 Roadmap
- [ ] Implement full MPI grid migration (octree re-partitioning).
- [ ] Port MHD (Magnetohydrodynamics) modules.
- [ ] Implement full Radiative Transfer (RT) support.

---
*Developed as part of the RAMSES-2025 Migration Task.*
