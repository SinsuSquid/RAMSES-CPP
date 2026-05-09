cmake_minimum_required(VERSION 3.15)
project(RAMSES-CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configuration Options
set(RAMSES_NDIM 3 CACHE STRING "Dimensions (1, 2, or 3)")
set(RAMSES_PRECISION 8 CACHE STRING "Real precision (4 or 8)")
set(RAMSES_NENER 0 CACHE STRING "Number of non-thermal energies")
set(RAMSES_NPSCAL 0 CACHE STRING "Number of passive scalars")
set(RAMSES_NMETALS 0 CACHE STRING "Number of metals (passive scalars)")
set(RAMSES_NGROUPS 0 CACHE STRING "Number of RT groups")
set(RAMSES_NIONS 0 CACHE STRING "Number of RT ions")
set(RAMSES_NX 1 CACHE STRING "Number of coarse cells in X")
set(RAMSES_NY 1 CACHE STRING "Number of coarse cells in Y")
set(RAMSES_NZ 1 CACHE STRING "Number of coarse cells in Z")
set(RAMSES_LONGINT OFF CACHE BOOL "Use 64-bit integers for grid IDs")
set(RAMSES_USE_MPI OFF CACHE BOOL "Enable MPI support")
set(RAMSES_USE_MHD OFF CACHE BOOL "Enable MHD solver")
set(RAMSES_USE_RT OFF CACHE BOOL "Enable RT solver")

# Source files for the library
set(LIB_SRCS
    src/AmrGrid.cpp
    src/Config.cpp
    src/Constants.cpp
    src/CoolingSolver.cpp
    src/HydroSolver.cpp
    src/RhdSolver.cpp
    src/Initializer.cpp
    src/RamsesReader.cpp
    src/RamsesWriter.cpp
    src/Simulation.cpp
    src/TreeUpdater.cpp
    src/Hilbert.cpp
    src/RiemannSolver.cpp
    src/RelativisticRiemannSolver.cpp
    src/SlopeLimiter.cpp
    src/Muscl.cpp
    src/PoissonSolver.cpp
    src/MhdSolver.cpp
    src/RtSolver.cpp
    src/RtChemistry.cpp
    src/MpiManager.cpp
    src/Parameters.cpp
    src/ParticleSolver.cpp
    src/ParticleSystem.cpp
    src/Cosmology.cpp
    src/LoadBalancer.cpp
    src/SolverFactory.cpp
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
    if(DEFINED RAMSES_${param})
        list(APPEND PARAM_DEFS "RAMSES_${param}=${RAMSES_${param}}")
    endif()
endforeach()

include_directories(include)

# Build a separate library and executable for each dimension
foreach(d 1 2 3)
    set(LIB_NAME ramses_lib_${d}d)
    add_library(${LIB_NAME} STATIC ${LIB_SRCS} ${CUSTOM_PATCH_SRCS})
    target_compile_definitions(${LIB_NAME} PUBLIC NDIM=${d} ${PARAM_DEFS})
    
    if(RAMSES_USE_MPI)
        find_package(MPI REQUIRED)
        target_link_libraries(${LIB_NAME} PUBLIC MPI::MPI_CXX)
        target_compile_definitions(${LIB_NAME} PUBLIC RAMSES_USE_MPI)
    endif()
    if(RAMSES_USE_MHD)
        target_compile_definitions(${LIB_NAME} PUBLIC MHD)
    endif()
    if(RAMSES_USE_RT)
        target_compile_definitions(${LIB_NAME} PUBLIC RT)
    endif()

    add_executable(ramses_${d}d src/main.cpp)
    target_compile_definitions(ramses_${d}d PRIVATE NDIM=${d})
    target_link_libraries(ramses_${d}d PRIVATE ${LIB_NAME})
endforeach()

# Reference verification tool
add_executable(verify_ref src/verify_ref.cpp)
target_compile_definitions(verify_ref PRIVATE NDIM=3)
target_link_libraries(verify_ref PRIVATE ramses_lib_3d)
