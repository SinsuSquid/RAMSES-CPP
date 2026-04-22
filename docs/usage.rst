Usage
=====

Running a Simulation
--------------------

RAMSES-CPP relies on standard Fortran Namelist (``.nml``) files to configure the simulation.

To run a simulation, specify the path to your namelist file as the first argument to the executable:

.. code-block:: bash

    ./ramses_main path/to/namelist.nml

For example, to run the provided 3D Sedov blast test:

.. code-block:: bash

    ./ramses_main ../namelist/sedov3d.nml

Parallel Execution
------------------

If you compiled RAMSES-CPP with MPI support, you can execute it across multiple processors using `mpirun` or `mpiexec`:

.. code-block:: bash

    mpirun -np 4 ./ramses_main ../namelist/sedov3d.nml

Output and Visualization
------------------------

Simulations output their state to directories named ``output_XXXXX/`` (where ``XXXXX`` is the output step number). 

Because RAMSES-CPP maintains strict binary compatibility with legacy RAMSES data formats, you can use existing Python visualization tools, such as the included ``visu_ramses.py``, to parse and plot these results seamlessly.
