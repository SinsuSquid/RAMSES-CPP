#!/bin/bash
#
# RAMSES-CPP Build Configuration Helper
# Adapted from legacy RAMSES build scripts for C++ version
#
# Usage:
#   ./build_config.sh [OPTIONS]
#
# Examples:
#   ./build_config.sh               # Default: MPI=OFF, MHD=ON, RT=ON
#   ./build_config.sh --mpi         # With MPI support
#   ./build_config.sh --no-mhd      # Without MHD
#   ./build_config.sh --mpi --no-rt # With MPI, without RT
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="${PROJECT_DIR}/build"

# Default configuration
USE_MPI="OFF"
USE_MHD="ON"
USE_RT="ON"
BUILD_TYPE="Release"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --mpi)
            USE_MPI="ON"
            shift
            ;;
        --no-mpi)
            USE_MPI="OFF"
            shift
            ;;
        --mhd)
            USE_MHD="ON"
            shift
            ;;
        --no-mhd)
            USE_MHD="OFF"
            shift
            ;;
        --rt)
            USE_RT="ON"
            shift
            ;;
        --no-rt)
            USE_RT="OFF"
            shift
            ;;
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --help)
            echo "RAMSES-CPP Build Configuration Helper"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --mpi, --no-mpi        Enable/disable MPI (default: OFF)"
            echo "  --mhd, --no-mhd        Enable/disable MHD (default: ON)"
            echo "  --rt, --no-rt          Enable/disable RT (default: ON)"
            echo "  --debug                Use Debug build type (default: Release)"
            echo "  --help                 Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure
echo "Configuring RAMSES-CPP..."
echo "  MPI:        $USE_MPI"
echo "  MHD:        $USE_MHD"
echo "  RT:         $USE_RT"
echo "  Build Type: $BUILD_TYPE"
echo ""

cmake .. \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DRAMSES_USE_MPI="$USE_MPI" \
    -DRAMSES_USE_MHD="$USE_MHD" \
    -DRAMSES_USE_RT="$USE_RT"

echo ""
echo "✅ Configuration complete! Run 'make -j1' to build."
