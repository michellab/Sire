#!/usr/bin/env bash

# Build script for Sire conda-forge package.

set -ex

# Make the build directories.
mkdir -p build/corelib
mkdir -p build/wrapper

# Build and install the core library.
cd build/corelib
cmake ${CMAKE_ARGS} \
    -D ANACONDA_BASE=${PREFIX} \
    -D BUILD_NCORES=$CPU_COUNT \
    -D SIRE_DISABLE_AVX=ON \
    -D SIRE_DISABLE_AVX512F=ON \
    ../../corelib
cmake --build . --target install -- VERBOSE=1 -j$CPU_COUNT

# Build and install the Python wrappers.
cd ../wrapper
cmake ${CMAKE_ARGS} \
    -D ANACONDA_BASE=${PREFIX} \
    -D BUILD_NCORES=$CPU_COUNT \
    ../../wrapper
cmake --build . --target install -- VERBOSE=1 -j$CPU_COUNT

# Remove the redundant files in the pkgs directory.
rm -r ${PREFIX}/pkgs/sire-*
