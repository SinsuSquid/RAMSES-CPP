<img src="logo/Sia-chan_logo.png" style="width: 900px; height: auto;">

Sia(_See-ah_)-chan, our mascot.

# RAMSES::CPP

A high-performance, modern C++17 port of the [**RAMSES-2025**](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Adaptive Mesh Refinement (AMR) code. 

RAMSES-CPP provides a functional, physically consistent alternative to the original Fortran implementation while maintaining strict binary compatibility for snapshots and restart files.

## 🚀 Status: Active Development
RAMSES-CPP is currently in the active development and verification phase. We are working towards:
- **Full Multi-Dimensional Hydrodynamics** (1D/2D/3D).
- **Magnetohydrodynamics (MHD)** stabilization in 2D/3D.
- **Gas Cooling and Heating** via an analytic ISM model.
- **Radiative Transfer (RT)** using the M1 closure scheme.
- **Self-Gravity** via a Multigrid Poisson solver.
- **Strict Binary Parity** with legacy Fortran snapshots (Under Verification).

For a detailed log of the migration progress, see [PORTING_HISTORY.md](PORTING_HISTORY.md).

---

## 📚 Documentation
Comprehensive documentation is available in the [docs/](docs/) directory:
- [**Installation**](docs/installation.md): How to build and configure the project.
- [**Usage**](docs/usage.md): Running simulations and using verification tools.
- [**Architecture**](docs/architecture.md): Overview of the C++ class structure and data flow.
- [**Physics Modules**](docs/index.md): Details on MHD, Gravity, and Hydro solvers.

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

2. **MHD Build:**
To enable MHD modules, use the `RAMSES_USE_MHD` flag:
```bash
cmake .. -DRAMSES_USE_MHD=ON
```

3. **Multi-Dimensional Build:**
To build for specific dimensions, use the `RAMSES_NDIM` flag:
```bash
cmake .. -DRAMSES_NDIM=1  # For 1D
cmake .. -DRAMSES_NDIM=2  # For 2D
```

---

## 🧪 Testing

The repository includes a comprehensive automated test suite to ensure physical consistency and binary parity with the legacy code. We are currently prioritizing the verification of:
- `hydro/sod-tube`: Classic shock tube benchmark.
- `hydro/sedov3d`: 3D blast wave benchmark.
- `mhd/orszag-tang`: 2D MHD vortex evolution.

To run the hydro suite:
```bash
cd tests
./run_test_suite.sh -t hydro
```

---

## 🏃 Running a Simulation

Use the standard RAMSES namelist format. The executable is located in your build directory (named `ramses_Nd` where N is the dimension):

```bash
# Example Sedov 3D test
./ramses_3d ../namelist/sedov3d.nml
```

Results are saved to `output_XXXXX/`, which is fully compatible with legacy RAMSES visualization tools.

---

## 🗺 Roadmap
- [x] Implement Hilbert-based domain decomposition and full state grid migration.
- [x] Implement dynamic MPI load balancing during simulation.
- [x] Port MHD (Magnetohydrodynamics) modules.
- [x] Implement gas cooling and heating (ISM model).
- [x] Implement basic Radiative Transfer (RT) support (M1 scheme).
- [x] Extend RT module with ionization and gas-coupling source terms.

## 🛡️ Verification and Parity

We ensure high-precision parity with legacy RAMSES using a dedicated C++ utility:

```bash
cd build
# Compare local snapshot with reference benchmark
./verify_ref ../tests/hydro/sod-tube/output_00001/amr_00001.out00001 output_00001/amr_00001.out00001
```

---
*Developed as part of the RAMSES-2025 Migration Task.*
