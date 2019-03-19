#!/usr/bin/env bash

# Create the archive name.
ARCHIVE=sire_conda_linux_latest.tar.bz2
if [ "$(uname)" == "Darwin" ]; then
    ARCHIVE=sire_conda_osx_latest.tar.bz2
fi

# Create a directory to store all of the necessary files for the Conda package.
mkdir sire_conda

# Create the required sub-directories.
mkdir -p sire_conda/lib/python3.7/site-packages

# Now copy the files into place.
cp -r $HOME/sire.app/pkgs/sire*/bin sire_conda
cp -r $HOME/sire.app/pkgs/sire*/include sire_conda
cp -r $HOME/sire.app/pkgs/sire*/lib sire_conda
cp -r $HOME/sire.app/pkgs/sire*/bundled/lib sire_conda
cp -r $HOME/sire.app/pkgs/sire*/share sire_conda

# Now copy across the Sire Python package.
cp -r $HOME/sire.app/lib/python3.7/site-packages/Sire sire_conda/lib/python3.7/site-packages

# Now compress the Conda package directory.
tar cjf $ARCHIVE sire_conda

# Remove the redundant package directory.
rm -rf sire_conda
