cmake_minimum_required(VERSION 3.15)
project(ramses-cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
set(RAMSES_NDIM 3 CACHE STRING "Number of dimensions")
set(RAMSES_PRECISION 8 CACHE STRING "Floating point precision")
set(RAMSES_NENER 0 CACHE STRING "Number of non-thermal energy variables")
option(RAMSES_LONGINT "Use 64-bit integers for particle IDs" ON)
option(RAMSES_USE_MPI "Enable MPI support" OFF)
option(RAMSES_USE_MHD "Enable MHD support" OFF)
option(RAMSES_USE_RT "Enable Radiative Transfer support" OFF)

# Add definitions
add_compile_definitions(NDIM=${RAMSES_NDIM})
add_compile_definitions(NENER=${RAMSES_NENER})
add_compile_definitions(RAMSES_PRECISION=${RAMSES_PRECISION})
if(RAMSES_LONGINT)
    add_compile_definitions(LONGINT)
endif()

if(RAMSES_USE_MHD)
    add_compile_definitions(MHD)
endif()

if(RAMSES_USE_RT)
    add_compile_definitions(RT)
endif()

if(RAMSES_USE_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    add_compile_definitions(RAMSES_USE_MPI)
endif()

# Include directories
include_directories(include)

# Source files
set(RAMSES_SOURCES
    src/Constants.cpp
    src/Parameters.cpp
    src/AmrGrid.cpp
    src/Config.cpp
    src/Hilbert.cpp
    src/TreeUpdater.cpp
    src/RamsesReader.cpp
    src/HydroSolver.cpp
    src/RiemannSolver.cpp
    src/SlopeLimiter.cpp
    src/Muscl.cpp
    src/Simulation.cpp
    src/MpiManager.cpp
    src/PoissonSolver.cpp
    src/Initializer.cpp
    src/RamsesWriter.cpp
    src/ParticleSystem.cpp
    src/ParticleSolver.cpp
    src/LoadBalancer.cpp
    src/MhdSolver.cpp
    src/CoolingSolver.cpp
    src/RtSolver.cpp
    src/RtChemistry.cpp
)

# Main library
add_library(ramses_lib STATIC ${RAMSES_SOURCES})
if(RAMSES_USE_MPI)
    target_link_libraries(ramses_lib PUBLIC MPI::MPI_CXX)
endif()

# Main executable
add_executable(ramses_${RAMSES_NDIM}d src/main.cpp)
target_link_libraries(ramses_${RAMSES_NDIM}d ramses_lib)

# Verification tool
add_executable(verify_ref src/verify_ref.cpp)
target_link_libraries(verify_ref ramses_lib)
