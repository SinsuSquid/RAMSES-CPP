cmake_minimum_required(VERSION 3.15)
project(ramses-cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
set(RAMSES_NDIM 3 CACHE STRING "Number of dimensions")
set(RAMSES_PRECISION 8 CACHE STRING "Floating point precision")
option(RAMSES_LONGINT "Use 64-bit integers for particle IDs" ON)

# Add definitions based on options
add_compile_definitions(NDIM=${RAMSES_NDIM})
add_compile_definitions(RAMSES_PRECISION=${RAMSES_PRECISION})
if(RAMSES_LONGINT)
    add_compile_definitions(LONGINT)
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
)

# Main library
add_library(ramses_lib STATIC ${RAMSES_SOURCES})

# Main executable (for now a placeholder)
add_executable(ramses_main src/main.cpp)
target_link_libraries(ramses_main ramses_lib)

# Testing (placeholder for GTest)
enable_testing()
# add_subdirectory(cpp/tests)
