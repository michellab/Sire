#!/usr/bin/env bash

# Set the source directory on macOS.
SRC_DIR=.

# Linux runs in a docker container from $HOME.
if [ ! -d $SRC_DIR/docker ]; then
    SRC_DIR=$HOME/Sire
fi

# Get the tag associated with the latest commit.
TAG=$(git --git-dir=$SRC_DIR/.git tag --contains)

# If the tag is empty, then simply exit
if [ -e $TAG ]; then
    exit 0
fi

# Convert . to _ in the tag.
TAG=$(echo $TAG | tr . _)

# Generate the name of the tagged binary release.
if [ "$(uname)" == "Linux" ]; then
    BINARY=sire_${TAG}_linux.run
else
    BINARY=sire_${TAG}_osx.run
fi

# Copy the latest devel binary to the release version.
cp $HOME/sire_devel_latest_*.run $HOME/$BINARY

# Deploy the binary.
$HOME/sire.app/bin/python $SRC_DIR/docker/sire-deploy-devel/deploy.py $HOME/$BINARY
