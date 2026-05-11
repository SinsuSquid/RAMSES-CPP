# RAMSES-CPP Project Instructions

This file defines the overarching conventions, architecture, and workflows for the RAMSES-CPP repository. All AI assistants and contributors must adhere to these guidelines.

## 1. Core Tech Stack & Build System
- **Language Standard:** C++17.
- **Build System:** CMake (>= 3.15).
- **Compile-Time Flags:** Ensure appropriate CMake options are used for dimensions (`RAMSES_NDIM`), precision (`RAMSES_PRECISION`), integer sizing (`RAMSES_LONGINT`), and physics modules (`RAMSES_USE_MPI`, `RAMSES_USE_MHD`, `RAMSES_USE_RT`).
- **MPI:** MPI support is optional. Default to sequential builds unless specifically testing multi-process behavior.

## 2. Code Style & Porting Patterns
- **Legacy Code Location:** The legacy Fortran RAMSES-2025 source code is located in the `./legacy` directory. Refer to this directory when porting logic or verifying behavior.
- **Modern C++ Idioms:** The project is a C++17 port of the Fortran RAMSES-2025 code. Translate Fortran idioms into safe, modern C++ patterns. Avoid C-style casts, unmanaged raw arrays (prefer `std::vector` or structured buffer abstractions), and magic numbers.
- **Architecture:** Maintain the modular solver class structure (e.g., `HydroSolver`, `MhdSolver`, `RtSolver`) as documented in `docs/architecture.md`.
- **Type Safety:** Do not bypass the type system or disable warnings. Use explicit instantiation and strong typing to maintain structural integrity.

## 3. Testing & Binary Parity
- **Test File Constraints:** You **MUST NOT** edit test-related data files located in `./tests/`, including namelists (`*.nml`), reference data (`*-ref.dat`), and plotting scripts (`plot-*.py`). However, you **ARE ALLOWED** to modify shell scripts (`*.sh`) in the `tests/` directory to adapt the testing infrastructure.
- **Strict Binary Compatibility:** Maintaining binary parity with the legacy RAMSES-2025 Fortran snapshots is a primary project goal.
- **Verification Workflow:** After making changes to physics or solver modules, verify snapshot parity using the included tool:
  ```bash
  cd build && ./verify_ref <reference_out_file> <local_out_file>
  ```
- **Automated Tests & Execution:** Always use a watchdog or timeout (e.g., `timeout 5m`) when running simulations or test suites to prevent stalled processes (due to numerical instability or tiny timesteps) from running indefinitely.
  ```bash
  export PYTHONPATH=${PYTHONPATH}:$(pwd)/tests/visu
  cd tests && timeout 10m ./run_test_suite.sh -t <suite_name>
  ```

## 4. Documentation & Version Control Workflows
- **Mandatory Updates:** Whenever there is an update or change to the codebase, you **MUST** ensure that a record is kept and corresponding updates are made in the following locations:
  - **`docs/` Directory:** Update relevant markdown guides (e.g., `architecture.md`, `usage.md`, or module-specific files) if system behavior changes.
  - **`README.md`:** Update the main README if there are changes to the project status, build instructions, or high-level overview.
  - **`PORTING_HISTORY.md`:** Log all significant migration steps, bug fixes related to parity, and newly implemented physics modules.
- **Git Operations:** After making code changes and updating the aforementioned records (`docs/`, `README.md`, `PORTING_HISTORY.md`), you are explicitly **allowed and expected to `git commit` and `git push`** the changes.
