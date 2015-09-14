#!/bin/bash

###
### This is a script that can be used to completely compile and
### install Sire from the current directory. This script will
### download and install miniconda, install openmm into that
### and will then configure and compile Sire in build/corelib
### and build/python
###

if [ -z "$INSTALL_SIRE_DIR" ]; then
    # Ask the user where they would like to install sire. By default
    # we will aim for $HOME/sire.app
    echo -n "Where would you like to install Sire? [$HOME/sire.app]: "
    read INSTALL_DIR

    if [ ! ${INSTALL_DIR} ]; then
        INSTALL_DIR=$HOME/sire.app
    else
        # Use eval so that we can expand variables such as $HOME
        INSTALL_DIR=`eval echo ${INSTALL_DIR}`
    fi
else
    INSTALL_DIR=`eval echo ${INSTALL_SIRE_DIR}`
fi

echo "Installing into directory '${INSTALL_DIR}'"

# Check whether or not miniconda has been downloaded and 
# installed into this directory
if [ -e "${INSTALL_DIR}/bin/python" ]; then
    # the miniconda distribution already exists.
    # We can now jump straight to the python install script
    ${INSTALL_DIR}/bin/python build/build_sire.py
    exit 0
fi

# Work out whether we are on OS X, Linux or Windows, and
# then download the appropriate miniconda distribution
if [ "$(uname)" == "Darwin" ]; then
    # This is running on a Mac
    PLATFORM="OSX"
    MINICONDA="http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # This is running on Linux
    PLATFORM="Linux"
    MINICONDA="http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    # This is running on Windows NT
    echo "Compilation on windows is not supported."
    exit -1
fi

# Download miniconda if it is not already in build/downloads
if [ ! -e "build/miniconda.sh" ]; then
    echo "Downloading miniconda from ${MINICONDA}..."
    wget ${MINICONDA} -O build/miniconda.sh
fi

# Now unpack miniconda and install it into the requested directory
echo "Unpacking miniconda into ${INSTALL_DIR}..."
bash build/miniconda.sh -b -p ${INSTALL_DIR}

# Now run the python install script
if [ -e "${INSTALL_DIR}/bin/python" ]; then
    ${INSTALL_DIR}/bin/python build/build_sire.py
    exit $?
else
    echo "** FATAL **"
    echo "** Cannot find ${INSTALL_DIR}/bin/python **"
    echo "** Something went wrong with the miniconda install! **"
    exit -1
fi
