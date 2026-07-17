cmake_minimum_required(VERSION 3.15)
project(RAMSES-CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build for performance (9.2x speedup over Debug)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.12.0
)
FetchContent_MakeAvailable(spdlog)

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

# Build library and executable for specified dimension only (from config.txt or default)
set(d ${NDIM})
set(LIB_NAME ramses_lib_${d}d)

# Initialize the library target
add_library(${LIB_NAME} STATIC)

# Custom Patches Discovery
file(GLOB CUSTOM_PATCH_SRCS "patch/*.cpp")
if(CUSTOM_PATCH_SRCS)
    message(STATUS "[RAMSES-CPP] Found custom C++ patches: ${CUSTOM_PATCH_SRCS}")
    add_definitions(-DRAMSES_USE_CUSTOM_PATCHES)
    target_sources(${LIB_NAME} PRIVATE ${CUSTOM_PATCH_SRCS})
endif()

# Traverse subdirectories to add source files
add_subdirectory(src)

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

target_compile_definitions(${LIB_NAME} PUBLIC NDIM=${d} NPRE=${NPRE} ${PARAM_DEFS})
target_link_libraries(${LIB_NAME} PUBLIC spdlog::spdlog)

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

# Dynamic git information gathering
string(TIMESTAMP BUILD_DATE "%m/%d/%y-%H:%M:%S")

# Build command line mock or approximate command
set(BUILD_CMD "cmake .. -DNDIM=${d} -DMPI=${MPI}")

# Retrieve Git hash, branch, and remote URL
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
execute_process(
    COMMAND git config --get remote.origin.url
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_REPO
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(NOT GIT_HASH)
    set(GIT_HASH "unknown")
endif()
if(NOT GIT_BRANCH)
    set(GIT_BRANCH "unknown")
endif()
if(NOT GIT_REPO)
    set(GIT_REPO "unknown")
endif()

target_compile_definitions(ramses_${d}d PRIVATE
    RAMSES_BUILD_DATE="${BUILD_DATE}"
    RAMSES_BUILD_COMMAND="${BUILD_CMD}"
    RAMSES_PATCH_DIR=""
    RAMSES_GIT_REPO="${GIT_REPO}"
    RAMSES_GIT_BRANCH="${GIT_BRANCH}"
    RAMSES_GIT_HASH="${GIT_HASH}"
)

target_link_libraries(ramses_${d}d PRIVATE ${LIB_NAME})


