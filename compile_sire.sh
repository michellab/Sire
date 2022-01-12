#!/bin/bash

###
### This is a script that can be used to completely compile and
### install Sire from the current directory. This script will
### download and install miniconda, install openmm into that
### and will then configure and compile Sire in build/corelib
### and build/python
###

build_and_install_sire()
{
# Now set the MACOSX_DEPLOYMENT_TARGET to make sure
# that we can work with Mavericks or above (needed by Qt5)
if [[ "$OSTYPE" == "darwin"* ]]; then
    export MACOSX_DEPLOYMENT_TARGET="10.9"
    export SDKROOT=`xcrun --show-sdk-path`
fi

# Now run the python install script
[ -z "$CONDA_PYTHON" ] && [ -e "${INSTALL_DIR}/bin/python" ] && CONDA_PYTHON="${INSTALL_DIR}/bin/python"
[ -z "$CONDA_PYTHON" ] && [ -e "${INSTALL_DIR}/python" ] && CONDA_PYTHON="${INSTALL_DIR}/python"
if [ ! -z "$CONDA_PYTHON" ]; then
    CONDA_BINDIR="`dirname "$CONDA_PYTHON"`"
    echo "** Running the conda activate script... **"
    echo "** . \"$CONDA_BINDIR/activate\""
    . "$CONDA_BINDIR/activate" ""
    echo "** Running the Python install script... **"
    echo "** \"$CONDA_PYTHON\" build/build_sire.py **"
    "$CONDA_PYTHON" build/build_sire.py
    err=$?
    conda deactivate
    exit $err
else
    echo "** FATAL **"
    echo "** Cannot find conda python **"
    echo "** Something went wrong with the miniconda install! **"
    echo "** Remove ${INSTALL_DIR}, then run compile_sire.sh --clean, then try again **"
    exit -1
fi
}

# Process arguments
key="$1"

case $key in
    -h|--help)
    echo "compile_sire.sh -h/--help shows help."
    echo "compile_sire.sh --install /path/to/sire.app will install Sire in /path/to/sire.app"
    echo "compile_sire.sh --clean completely cleans the build directory."
    exit 0
    ;;
    --install)
    INSTALL_SIRE_DIR="$2"
    echo "Installing Sire into ${INSTALL_SIRE_DIR}"
    ;;
    --clean)
    echo "Completely cleaning the build directories..."
    echo "rm -rf build/miniconda.sh build/corelib/* build/wrapper/*"
    rm -rf build/miniconda.sh build/corelib/* build/wrapper/*
    echo "...all clean"
    exit 0
    ;;
esac

# Set the version of miniconda to use. Choose "latest" for the latest
# miniconda, or set a specific version here
MINICONDA_VERSION="4.10.1"
PYTHON_VERSION="py38"
#MINICONDA_VERSION="latest"

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
if [ -e "${INSTALL_DIR}/bin/python" ] || [ -e "${INSTALL_DIR}/python.exe" ]; then
    build_and_install_sire
    exit 0
fi

# Now work out if we are 32bit or 64bit...
LONG_BIT=$(getconf LONG_BIT)
OS_BIT=$(uname -m)

# (note that BIT_TYPE is the string used by miniconda to
#  denote a 32bit or 64bit operating system)
BIT_TYPE="x86_64"

if [ ${LONG_BIT} == "64" ]; then
    echo "64 bit processor"
elif [ ${LONG_BIT} == "32" ]; then
    echo "32 bit processor"
else
    echo "Unknown number of bits. Assuming 64 bit"
fi

if [ ${OS_BIT} == "x86_64" ]; then
    echo "x86_64 operating system"
elif [ ${OS_BIT} == "ppc64le" ]; then
    echo "ppc64le operating system"
    BIT_TYPE=ppc64le
elif [ ${OS_BIT} == "arm64" ]; then
    echo "arm64 operating system"
    BIT_TYPE="arm64"
elif [ ${LONG_BIT} == "32" ]; then
    # must be 32 bit operating system on 32 bit processor
    echo "i386 operating system"
    BIT_TYPE="x86"
elif [ ${OS_BIT} == "i386" ] || [ ${OS_BIT} == "i486" ] || [ ${OS_BIT} == "i586" ] || [ ${OS_BIT} == "i686" ]; then
    echo "i386 operating system"
    BIT_TYPE="x86"
else
    echo "Assuming a 64 bit operating system."
    echo "If you want to change this, set BIT_TYPE to x86 here"
    #BIT_TYPE="x86"
fi

# Work out whether we are on OS X, Linux or Windows, and
# then download the appropriate miniconda distribution
if [ "$(uname)" == "Darwin" ]; then
    # This is running on a Mac
    PLATFORM="OSX"
    MINICONDA="https://repo.continuum.io/miniconda/Miniconda3-${PYTHON_VERSION}_${MINICONDA_VERSION}-MacOSX-${BIT_TYPE}.sh"
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # This is running on Linux
    PLATFORM="Linux"
    MINICONDA="https://repo.continuum.io/miniconda/Miniconda3-${PYTHON_VERSION}_${MINICONDA_VERSION}-Linux-${BIT_TYPE}.sh"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    # This is running on Windows NT
    echo "Compilation on windows is not supported."
    exit -1
elif [ "$(expr substr $(uname -s) 1 9)" == "CYGWIN_NT" ]; then
    # This is running on windows under cygwin
    echo "Running an install under cygwin on windows"
    PLATFORM="Windows"
    SUBPLATFORM="Cygwin"
    MINICONDA="https://repo.continuum.io/miniconda/Miniconda3-${PYTHON_VERSION}_${MINICONDA_VERSION}-Windows-${BIT_TYPE}.exe"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    # This is running on windows under cygwin
    echo "Running an install under MSYS2 on windows"
    PLATFORM="Windows"
    SUBPLATFORM="MSYS2"
    MINICONDA="https://repo.continuum.io/miniconda/Miniconda3-${PYTHON_VERSION}_${MINICONDA_VERSION}-Windows-${BIT_TYPE}.exe"
else
    # Cannot identify the platform. Tell the user to download
    # miniconda directly
    echo "Cannot identify your operating system / platform."
    echo "Please download and install miniconda manually, and then "
    echo "use the Python from that miniconda to run the build/build_sire.py script."
    PLATFORM="$(uname -s)"
    echo $PLATFORM
    exit -1
fi

if [ -e ${INSTALL_DIR} ]; then
    echo "Install directory already exists. Assuming that miniconda is already installed here."
else
    if [ ${PLATFORM} == "Windows" ]; then
        # Download miniconda.exe if it is not already in build/downloads
        if [ ! -e "build/miniconda.exe" ]; then
            echo "** Downloading miniconda from ${MINICONDA}... **"
            echo "curl -L ${MINICONDA} -o build/miniconda.exe"
            curl -L ${MINICONDA} -o build/miniconda.exe
            chmod a+x build/miniconda.exe
        fi

        # Now unpack miniconda and install it into the requested directory
        echo "Running the miniconda installation. Make sure you install miniconda just for yourself."
        echo "Also, ensure that you install miniconda into the directory 'C:\msys2\${INSTALL_DIR}'"
        echo "Also note that you should't select the option to 'add anaconda to the PATH' or to"
        echo "register anaconda as the default python"
        ./build/miniconda.exe
    else
        # Download miniconda if it is not already in build/downloads
        if [ ! -e "build/miniconda.sh" ]; then
            echo "** Downloading miniconda from ${MINICONDA}... **"
            echo "** wget ${MINICONDA} -O build/miniconda.sh 2>/dev/null || curl -L ${MINICONDA} -o build/miniconda.sh  **"
            wget ${MINICONDA} -O build/miniconda.sh 2>/dev/null || curl -L ${MINICONDA} -o build/miniconda.sh
        fi

        # Now unpack miniconda and install it into the requested directory
        echo "** Unpacking miniconda into ${INSTALL_DIR}... **"
        echo "** bash build/miniconda.sh -b -p ${INSTALL_DIR} **"
        bash build/miniconda.sh -b -p ${INSTALL_DIR}
    fi
fi

build_and_install_sire
