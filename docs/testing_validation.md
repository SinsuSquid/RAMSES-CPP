# Testing & Validation Protocol

This document describes the design of the differential test suite, the python verification pipeline, and the step-by-step procedure for validating binary parity.

---

## 🧪 Differential Testing Suite

To guarantee that the modern C++ engine matches the original Fortran physics, RAMSES-CPP uses a **differential testing suite**. 

* **The Core Mechanism:** For each physical test, the engine runs a simulation and compares the resulting output snapshots cell-by-cell against pre-compiled legacy references.
* **The Test Runner:** [run_test_suite.sh](file:///home/bgkang/Projects/RAMSES-CPP/tests/run_test_suite.sh) automatically configures CMake, builds the target dimensional executable (automatically detecting and leveraging `ninja` if installed on your system for faster incremental builds), runs the simulation, and executes the Python validation checks.

### Test Runner Options
```bash
# Run all tests (requires timeout to prevent hangs)
timeout 10m ./run_test_suite.sh

# Run only hydrodynamics tests
timeout 10m ./run_test_suite.sh -t hydro

# Run with MPI support on 4 processors
timeout 10m ./run_test_suite.sh -p 4 -t mhd

# Run with code coverage analysis enabled
timeout 10m ./run_test_suite.sh -s -t poisson
```

---

## 📁 Anatomy of a Test Case

Each test folder under `tests/<category>/<test_name>/` contains:
1. **`config.txt`:** Records compilation parameters needed for the test:
   `FLAGS: NDIM=1 PATCH= SOLVER=hydro`
2. **`<test_name>.nml`:** The standard RAMSES namelist containing run-time parameters (resolution, box length, fluid equations).
3. **`plot-<test_name>.py`:** The plotting and analysis script. It reads local data, compares it cell-by-cell to reference data, generates PDF plots, and prints `PASSED` or `FAILED`.
4. **`<test_name>-ref.dat`:** Pre-compiled reference snapshot values containing cell coords, density, velocity, and pressure.

---

## 🐍 Python Parity Validation (`visu_ramses.py`)

Snapshot data comparison is managed by the python module [visu_ramses.py](file:///home/bgkang/Projects/RAMSES-CPP/tests/visu/visu_ramses.py).

### How Comparison Works
1. The python script loads the C++ binary output files (using custom binary readers that mirror Fortran unformatted record structures).
2. It extracts physical arrays (density, velocity, pressure, magnetic fields, etc.) for all cells.
3. It compares these values cell-by-cell against the reference snapshot:
   $$\max \left( |U_{\text{local}} - U_{\text{ref}}| \right) < \text{tolerance}$$
   * The default tolerance is extremely strict (typically $3 \cdot 10^{-12}$ for double-precision runs).
4. If the difference is within the tolerance, it prints `PASSED` and the test runner reports `[ OK ]`. Otherwise, it prints `FAILED` and reports `[FAIL]`.

---

## ➕ How to Add a New Test

To add a new physics regression test:

1. **Create the Folder:** Create a sub-folder under the appropriate category, e.g., `tests/hydro/my_test/`.
2. **Add Namelist & Config:** Add your `my_test.nml` and a `config.txt` specifying the compilation flags.
3. **Generate Reference Data:**
   * Run the legacy Fortran RAMSES code for your namelist.
   * Extract cell coordinates and variable values to a text/binary file and save it as `my_test-ref.dat` in the folder.
4. **Write the Plotting Script:**
   * Create `plot-my_test.py`.
   * Load local output snapshots and reference snapshots.
   * Add comparison assertions (density, pressure, etc.) and output a plot.
   * Ensure it prints `PASSED` on success.
5. **Update Test List:** The test runner scans directories automatically. Running `./run_test_suite.sh` will discover and execute your new test.
