#!/usr/bin/env bash

# Create the archive name.
ARCHIVE=sire_conda_linux.tar.bz2
if [ "$(uname)" == "Darwin" ]; then
    ARCHIVE=sire_conda_osx.tar.bz2
fi

# Create a directory to store all of the necessary files for the Conda package.
mkdir sire_conda

# Create the required sub-directories.
mkdir sire_conda/bin
mkdir -p sire_conda/lib/python3.7/site-packages
mkdir sire_conda/pkgs

# Now copy the files into place.

# First the binary files.
for i in $(ls -l $HOME/sire.app/bin | grep sire_python | awk '{print $9}')
do
    cp $HOME/sire.app/bin/$i sire_conda/bin
done

# Next the Sire library files.
cp $HOME/sire.app/lib/libSire* sire_conda/lib
cp $HOME/sire.app/lib/libSquire* sire_conda/lib

# Add the bundled libcpuid files.
cp $HOME/sire.app/lib/libcpuid* sire_conda/lib

# Now copy across the Sire Python package.
cp -r $HOME/sire.app/lib/python3.7/site-packages/Sire sire_conda/lib/python3.7/site-packages

# Finally the Sire pkg files.
cp -r $HOME/sire.app/pkgs/sire* sire_conda/pkgs

# Now compress the Conda package directory.
tar -czf $ARCHIVE sire_conda

# Remove the redundant package directory.
rm -rf sire_conda
