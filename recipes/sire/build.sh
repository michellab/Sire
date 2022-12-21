#!/usr/bin/env bash

# Build script for Sire conda-forge package.

# Exit immediately on error.
set -e

# Build/install Sire.
python setup.py install --skip-deps

# Remove the build files to free up space.
rm -r build
