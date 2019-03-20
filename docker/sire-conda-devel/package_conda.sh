#!/usr/bin/env bash

# Create the archive name.
ARCHIVE=sire_conda_latest_linux.tar.bz2
if [ "$(uname)" == "Darwin" ] ; then
    ARCHIVE=sire_conda_latest_osx.tar.bz2
fi

# Create a directory to store all of the necessary files for the Conda package.
mkdir sire_conda

# Create the required sub-directories.
mkdir -p sire_conda/lib/python3.7/site-packages

# Copy the files into place.
cp -r $HOME/sire.app/pkgs/sire*/bin sire_conda
cp -r $HOME/sire.app/pkgs/sire*/include sire_conda
cp -r $HOME/sire.app/pkgs/sire*/lib sire_conda
cp -r $HOME/sire.app/pkgs/sire*/bundled/lib sire_conda
cp -r $HOME/sire.app/pkgs/sire*/share sire_conda

# Copy across the Sire Python package.
cp -r $HOME/sire.app/lib/python3.7/site-packages/Sire sire_conda/lib/python3.7/site-packages

# Compress the Conda package directory.
tar cjf $ARCHIVE sire_conda

# Remove the redundant package directory.
rm -rf sire_conda
