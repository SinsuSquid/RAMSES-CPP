---
layout: default
title: Usage
---

# Usage

RAMSES-CPP is designed to be plug-and-play with the existing ecosystem of the original Fortran RAMSES code.

## Running a Simulation

The simulation is configured using standard Fortran Namelist (`.nml`) files.

To run a simulation, execute `ramses_main` and pass the path to your namelist file as the first argument:

```bash
cd build
./ramses_main path/to/your/namelist.nml
```

### Example: 3D Sedov Blast

A standard Sedov 3D blast wave test is included in the `namelist` directory:

```bash
./ramses_main ../namelist/sedov3d.nml
```

## Parallel Execution (MPI)

If you compiled RAMSES-CPP with MPI support, you can execute it across multiple processors using `mpirun` or `mpiexec`:

```bash
mpirun -np 4 ./ramses_main ../namelist/sedov3d.nml
```

The code will automatically perform domain decomposition using a Hilbert curve and distribute the workload across the available ranks.

## Output and Visualization

Simulations output their state to directories named `output_XXXXX/` (where `XXXXX` is the output step number). 

Because RAMSES-CPP maintains **strict binary compatibility** with legacy RAMSES data formats, it generates the expected `amr`, `hydro`, and `info` files. 

You can use existing Python visualization tools to parse and plot these results seamlessly:

```bash
# Example using the included Python plotting script for a hydro test
python3 ../tests/hydro/implosion/plot-implosion.py
```

## Automated Testing

The repository includes a comprehensive test suite to verify the physics and AMR logic.

To run the tests:
```bash
cd tests
./run_test_suite.sh -t hydro  # Run hydrodynamics tests
./run_test_suite.sh -t mhd    # Run MHD tests (requires MHD build)
```

The test suite automatically compares the C++ output against reference data and reports any discrepancies.

## Verification Tool

The repository includes a `verify_ref` utility to compare C++ snapshots against reference Fortran results to ensure bit-perfect structural parity.

```bash
./verify_ref path/to/fortran_amr path/to/cpp_amr
```
