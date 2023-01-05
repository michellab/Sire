#!/usr/bin/env bash

# Build script for Sire conda-forge package.

# Exit immediately on error.
set -e

# Build/install Sire.
python setup.py install --skip-deps
