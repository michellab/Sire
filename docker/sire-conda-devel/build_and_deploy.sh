#!/usr/bin/env bash

# Get the Anaconda usernam and access token.
ANACONDA_TOKEN=${ANACONDA_TOKEN:-$1}

# Check the OS.
if [ "$(uname)" == "Darwin" ]; then
    OS=osx-64
else
    OS=linux-64
fi

# Set the bin directory.
BIN_DIR=$HOME/sire.app/bin

# Set the Conda build directory.
CONDA_DIR=$HOME/conda

# Move the to build directory.
cd $CONDA_DIR

# Build the Conda package.
$BIN_DIR/conda-build -c conda-forge .

# Upload the package to the michellab channel on anaconda.org
$BIN_DIR/anaconda -t $ANACONDA_TOKEN --user michellab $HOME/sire.app/conda-bld/$OS/sire-* --label dev
