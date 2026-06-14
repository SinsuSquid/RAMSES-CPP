cmake_minimum_required(VERSION 3.15)
project(RAMSES-CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build for performance (9.2x speedup over Debug)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Configuration Options
set(NDIM 3 CACHE STRING "Dimensions (1, 2, or 3)")
set(NPRE 8 CACHE STRING "Real precision (4 or 8)")
set(NENER 0 CACHE STRING "Number of non-thermal energies")
set(NPSCAL 0 CACHE STRING "Number of passive scalars")
set(NMETALS 0 CACHE STRING "Number of metals (passive scalars)")
set(NGROUPS 0 CACHE STRING "Number of RT groups")
set(NIONS 0 CACHE STRING "Number of RT ions")
set(NX 1 CACHE STRING "Number of coarse cells in X")
set(NY 1 CACHE STRING "Number of coarse cells in Y")
set(NZ 1 CACHE STRING "Number of coarse cells in Z")
set(LONGINT OFF CACHE BOOL "Use 64-bit integers for grid IDs")
set(MPI OFF CACHE BOOL "Enable MPI support")
set(SOLVER "hydro" CACHE STRING "Solver type (hydro, mhd, or rhd)")
set(RT OFF CACHE BOOL "Enable RT solver")
set(ATON OFF CACHE BOOL "Enable ATON GPU solver")
set(GRACKLE OFF CACHE BOOL "Enable Grackle cooling library")
set(USE_TURB OFF CACHE BOOL "Enable FFTW-based turbulence driving")
set(NVECTOR 32 CACHE STRING "Size of vector cache")

# Source files for the library
set(LIB_SRCS
    # Core
    src/core/AmrGrid.cpp
    src/core/Config.cpp
    src/core/Constants.cpp
    src/core/Parameters.cpp
    src/core/MpiManager.cpp
    src/core/Simulation.cpp
    src/core/TreeUpdater.cpp
    src/core/Hilbert.cpp
    src/core/LoadBalancer.cpp
    src/core/SolverFactory.cpp
    src/core/Cosmology.cpp
    
    # Solvers
    src/solvers/hydro/HydroSolver.cpp
    src/solvers/hydro/RiemannSolver.cpp
    src/solvers/hydro/SlopeLimiter.cpp
    src/solvers/mhd/MhdSolver.cpp
    src/solvers/rhd/RhdSolver.cpp
    src/solvers/rhd/RelativisticRiemannSolver.cpp
    src/solvers/rt/RtSolver.cpp
    src/solvers/rt/RtChemistry.cpp
    src/solvers/gravity/PoissonSolver.cpp
    src/solvers/physics/CoolingSolver.cpp
    src/solvers/physics/StarSolver.cpp
    src/solvers/physics/FeedbackSolver.cpp
    src/solvers/physics/SinkSolver.cpp
    src/solvers/physics/TurbulenceSolver.cpp
    src/solvers/physics/ClumpFinder.cpp

    # Particles
    src/particles/ParticleSolver.cpp
    src/particles/ParticleSystem.cpp

    # IO / Diagnostics
    src/io/RamsesReader.cpp
    src/io/RamsesWriter.cpp
    src/io/Hdf5Writer.cpp
    src/io/LightCone.cpp
    
    # Root Initialization
    src/Initializer.cpp
)

# Custom Patches Discovery
file(GLOB CUSTOM_PATCH_SRCS "patch/*.cpp")
if(CUSTOM_PATCH_SRCS)
    message(STATUS "[RAMSES-CPP] Found custom C++ patches: ${CUSTOM_PATCH_SRCS}")
    add_definitions(-DRAMSES_USE_CUSTOM_PATCHES)
endif()

# Global definitions from parameters
set(PARAM_DEFS "")
foreach(param NX NY NZ NGROUPS NIONS NMETALS NENER NPSCAL)
    if(DEFINED ${param})
        list(APPEND PARAM_DEFS "RAMSES_${param}=${${param}}")
    endif()
endforeach()

include_directories(
    include
    include/ramses
    include/ramses/core
    include/ramses/solvers/hydro
    include/ramses/solvers/mhd
    include/ramses/solvers/rhd
    include/ramses/solvers/rt
    include/ramses/solvers/gravity
    include/ramses/solvers/physics
    include/ramses/particles
    include/ramses/io
)

# Build library and executable for specified dimension only (from config.txt or default)
set(d ${NDIM})
set(LIB_NAME ramses_lib_${d}d)
add_library(${LIB_NAME} STATIC ${LIB_SRCS} ${CUSTOM_PATCH_SRCS})
target_compile_definitions(${LIB_NAME} PUBLIC NDIM=${d} NPRE=${NPRE} ${PARAM_DEFS})

if(MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(${LIB_NAME} PUBLIC MPI::MPI_CXX)
    target_compile_definitions(${LIB_NAME} PUBLIC RAMSES_USE_MPI)
endif()
if(SOLVER STREQUAL "mhd")
    target_compile_definitions(${LIB_NAME} PUBLIC MHD)
endif()
if(SOLVER STREQUAL "rhd")
    target_compile_definitions(${LIB_NAME} PUBLIC RHD)
endif()
if(RT)
    target_compile_definitions(${LIB_NAME} PUBLIC RT)
endif()
if(ATON)
    target_compile_definitions(${LIB_NAME} PUBLIC ATON)
endif()
if(GRACKLE)
    target_compile_definitions(${LIB_NAME} PUBLIC grackle)
endif()
if(USE_TURB)
    target_compile_definitions(${LIB_NAME} PUBLIC USE_TURB)
endif()

add_executable(ramses_${d}d src/main.cpp)
target_compile_definitions(ramses_${d}d PRIVATE NDIM=${d})
target_link_libraries(ramses_${d}d PRIVATE ${LIB_NAME})


