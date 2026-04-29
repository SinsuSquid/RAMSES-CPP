<img src="logo/Sia-chan_logo.png" style="width: 900px; height: auto;">

Sia(_See-ah_)-chan, our mascot.

# RAMSES::CPP

A high-performance, modern C++17 port of the [**RAMSES-2025**](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Adaptive Mesh Refinement (AMR) code. 

RAMSES-CPP provides a functional, physically consistent alternative to the original Fortran implementation while maintaining strict binary compatibility for snapshots and restart files.

## 🚀 Status: Production-Ready
The code now supports:
- **Full Multi-Dimensional Support:** Seamlessly build for 1D, 2D, or 3D using CMake flags.
- **3D Unsplit Godunov Hydrodynamics:** HLLC/LLF Riemann solvers with MUSCL-Hancock reconstruction.
- **Multigrid Poisson Solver:** Iterative Gauss-Seidel for self-gravity.
- **N-Body Dynamics:** Full Particle-Mesh (CIC) assignment and advection.
- **AMR Engine:** Recursive sub-cycling with dynamic tree traversal.
- **Ecosystem Parity:** Strict unformatted Fortran binary parity, compatible with legacy RAMSES tools (`visu_ramses.py`).

---

## 🛠 Building the Project

### Prerequisites:
- CMake (>= 3.15)
- C++17 compliant compiler (GCC 9+, Clang, etc.)
- MPI (Optional: detect automatically, defaults to sequential if absent)

### Build Instructions:

1. **Standard Build (3D):**
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

2. **Multi-Dimensional Build:**
To build for specific dimensions, use the `RAMSES_NDIM` flag:
```bash
cmake .. -DRAMSES_NDIM=1  # For 1D
cmake .. -DRAMSES_NDIM=2  # For 2D
```

---

## 🧪 Testing

The repository includes a comprehensive automated test suite to ensure physical consistency and binary parity with the legacy code.

```bash
# Run the full suite
cd tests
./run_test_suite.sh -t hydro
```

Verified tests include:
- `hydro/sod-tube`: Classic shock tube benchmark.
- `hydro/sod-tube-nener`: Advanced test with multiple passive energy variables.
- `hydro/barotrop`: Self-gravitating collapsing sphere.

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
