<img src="logo/Sia-chan.png" width="100%">

Sia(_See-ah_)-chan, our mascot.

# RAMSES::CPP

A high-performance, modern C++17 port of the [**RAMSES-2025**](https://github.com/ramses-organisation/ramses/releases/tag/2025.05) Adaptive Mesh Refinement (AMR) code. 

RAMSES-CPP provides a functional, physically consistent alternative to the original Fortran implementation while maintaining strict binary compatibility for snapshots and restart files.

## 🚀 Status: Phase 3 Complete (Scalar Verification)
RAMSES-CPP has successfully completed the scalar transport verification phase! We have achieved:
- [x] **Passive Scalar Parity:** Verified 100% correct advection of passive mass fractions, fully separated from gas thermodynamics.
- [x] **Memory Safety:** Leveraged Valgrind to eliminate critical out-of-bounds accesses in the AMR engine, ensuring production-grade stability.
- [x] **Native Root Grids:** Refactored the AMR data structure to natively support Level 1 root grids, ensuring bit-perfect binary parity and full box coverage in snapshots.
- [x] **Physical Instabilities:** Successfully seeded wake turbulence (Kelvin-Helmholtz) in 2D benchmarks using additive velocity jitter.
- [x] **Adaptive Mesh Refinement:** Upgraded the flagging logic to support de-refinement, enabling efficient mesh tracking of scalar gradients.
- [x] **Snapshot Bit-Perfection:** Achieved bit-perfect binary parity for AMR and Hydro snapshots across all dimensions, ensuring full compatibility with legacy visualization scripts.

We are now entering Phase 19: **RT Core Implementation (M1 Transport)**.

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

The repository includes a comprehensive automated test suite to ensure physical consistency and binary parity with the legacy code.
Ensure `./tests/visu/` is in your `PYTHONPATH` before running tests.

To run the hydro suite:
```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/tests/visu
cd tests
./run_test_suite.sh -t hydro
```

---

## 🏃 Running a Simulation

Use the standard RAMSES namelist format. The executable is located in your build directory (named `ramses_Nd` where N is the dimension):

```bash
# Example 1D Advection test
./ramses_1d ../namelist/advect1d.nml
```

Results are saved to `output_XXXXX/`, which is fully compatible with legacy RAMSES visualization tools.

---

## 🗺 Roadmap
- [x] Port MUSCL-Hancock 2nd-order hydro solver.
- [x] Implement gradient-based AMR refinement criteria.
- [x] Implement Hilbert-based domain decomposition and full state grid migration.
- [x] Implement dynamic MPI load balancing during simulation.
- [x] Port MHD (Magnetohydrodynamics) modules.
- [x] Implement gas cooling and heating (ISM model).
- [x] Implement basic Radiative Transfer (RT) support (M1 scheme).
- [ ] Extend RT module with ionization and gas-coupling source terms.
- [ ] Reach 100% test coverage for the `hydro` and `mhd` test suites.
