#!/bin/bash
NDIM=$1
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BASE_DIR}/build"
BIN_DIR="${BASE_DIR}/bin"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"
cmake .. -DRAMSES_NDIM=${NDIM} -DCMAKE_BUILD_TYPE=Debug
make ramses_main

mkdir -p "${BIN_DIR}"
cp ramses_main "${BIN_DIR}/test_exe_${NDIM}d"
