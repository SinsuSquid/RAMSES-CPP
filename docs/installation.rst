Installation
============

Prerequisites
-------------

To build RAMSES-CPP, you need the following:

- **CMake** (>= 3.15)
- **C++17 compliant compiler** (e.g., GCC 9+, Clang)
- **MPI** (Optional: automatically detected for parallel execution)

Building from Source
--------------------

1. Clone the repository and navigate to the project directory:

.. code-block:: bash

    git clone https://github.com/SinsuSquid/RAMSES_CPP.git
    cd RAMSES_CPP

2. Create a build directory and run CMake:

.. code-block:: bash

    mkdir build && cd build
    cmake ..

3. Compile the project:

.. code-block:: bash

    make -j$(nproc)

This will generate the `ramses_main` executable in the `build/` directory.

Parallel Support
----------------

If MPI is installed and detected by CMake, the build system will automatically define ``RAMSES_USE_MPI``, enabling the ``LoadBalancer`` for distributed execution.
