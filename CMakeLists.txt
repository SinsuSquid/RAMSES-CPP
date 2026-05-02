cmake_minimum_required(VERSION 3.15)
project(RAMSES-CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configuration Options
set(RAMSES_NDIM 3 CACHE STRING "Dimensions (1, 2, or 3)")
set(RAMSES_PRECISION 8 CACHE STRING "Real precision (4 or 8)")
set(RAMSES_LONGINT OFF CACHE BOOL "Use 64-bit integers for grid IDs")
set(RAMSES_USE_MPI OFF CACHE BOOL "Enable MPI support")
set(RAMSES_USE_MHD OFF CACHE BOOL "Enable MHD solver")
set(RAMSES_USE_RT OFF CACHE BOOL "Enable RT solver")

add_definitions(-DNDIM=${RAMSES_NDIM})

if(RAMSES_USE_MPI)
    find_package(MPI REQUIRED)
    add_definitions(-DRAMSES_USE_MPI)
endif()

if(RAMSES_USE_MHD)
    add_definitions(-DMHD)
endif()

if(RAMSES_USE_RT)
    add_definitions(-DRT)
endif()

include_directories(include)

# Core library
add_library(ramses_lib STATIC
    src/AmrGrid.cpp
    src/Config.cpp
    src/Constants.cpp
    src/HydroSolver.cpp
    src/Initializer.cpp
    src/RamsesReader.cpp
    src/RamsesWriter.cpp
    src/Simulation.cpp
    src/TreeUpdater.cpp
    src/Hilbert.cpp
    src/RiemannSolver.cpp
    src/SlopeLimiter.cpp
    src/Muscl.cpp
    src/PoissonSolver.cpp
    src/MpiManager.cpp
    src/Parameters.cpp
)

if(RAMSES_USE_MPI)
    target_link_libraries(ramses_lib PUBLIC MPI::MPI_CXX)
endif()

# Main executable
add_executable(ramses_${RAMSES_NDIM}d src/main.cpp)
target_link_libraries(ramses_${RAMSES_NDIM}d PRIVATE ramses_lib)

# Reference verification tool
add_executable(verify_ref src/verify_ref.cpp)
target_link_libraries(verify_ref PRIVATE ramses_lib)
