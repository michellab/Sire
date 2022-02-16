#!/usr/bin/env bash

set -ex

if [ "$(uname)" == "Linux" ]; then
    # Update the C compiler on Linux.
    sed -i "s#gcc#$CC#g" Makefile.fkcombu
else
    # Ignore implicit function declaration errors on macOS.
    sed -i.bak -e "s#Wall#Wno-implicit-function-declaration#g" Makefile.fkcombu
fi

# Build FKCOMBU.
make -f Makefile.fkcombu

# Copy the FKCOMBU executable to the bin directory.
mkdir ${PREFIX}/bin
cp -a ../fkcombu ${PREFIX}/bin
