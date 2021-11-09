#!/usr/bin/env bash

# Get the Anaconda Cloud access token.
ANACONDA_TOKEN=${ANACONDA_TOKEN:-$1}

# Check the OS.
if [ "$(uname)" == "Darwin" ]; then
    OS=osx-64
else
    OS=linux-64
fi

# Set the bin directory.
BIN_DIR=$HOME/sire.app/bin

# Set the source and Conda build directory on macOS.
SRC_DIR=$(pwd)
CONDA_DIR=$SRC_DIR/docker/sire-conda/recipe

# Linux runs in a docker container from $HOME.
if [ ! -d $CONDA_DIR ]; then
    SRC_DIR=$HOME/Sire
    CONDA_DIR=$HOME/Sire/docker/sire-conda/recipe
fi

# Move the to build directory.
cd $CONDA_DIR

# Set the default Conda label.
LABEL=dev

# Get the tag associated with the latest commit.
TAG=$(git --git-dir=$SRC_DIR/.git --work-tree=$SRC_DIR tag --contains)

# If the tag is not empty, then set the label to main.
if [ ! -e $TAG ]; then
    LABEL=main
fi

# Build the Conda package.
$BIN_DIR/conda-build -c conda-forge -c omnia .

# Upload the package to the michellab channel on Anaconda Cloud.

# Label release packages with main and dev so that dev is at least as new as main.
if [ "$LABEL" = "main" ]; then
    $BIN_DIR/anaconda \
        --token $ANACONDA_TOKEN upload \
        --user michellab \
        --label main \
        --label dev \
        $HOME/sire.app/conda-bld/$OS/sire-*
else
    $BIN_DIR/anaconda \
        --token $ANACONDA_TOKEN upload \
        --user michellab \
        --label dev \
        $HOME/sire.app/conda-bld/$OS/sire-*
fi
