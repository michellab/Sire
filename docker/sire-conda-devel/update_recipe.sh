#!/usr/bin/env bash

# Set the Conda build directory on macOS.
CONDA_DIR=./docker/sire-conda-devel/recipe

# Linux runs in a docker container from $HOME.
if [ ! -d $CONDA_DIR ]; then
    CONDA_DIR=$HOME/Sire/docker/sire-conda-devel/recipe
fi

# Set the Sire source directory.
if [ -d corelib ]; then
    SRC_DIR=.
else
    SRC_DIR=$HOME/Sire
fi

# Store the name of the recipe and template YAML files.
RECIPE=$CONDA_DIR/meta.yaml
TEMPLATE=$CONDA_DIR/template.yaml

# Overwite the recipe with the template file.
cp $TEMPLATE $RECIPE

ARCHIVE_LINUX=$HOME/sire_conda_latest_linux.tar.bz2
ARCHIVE_OSX=$HOME/sire_conda_latest_osx.tar.bz2

# List of python dependencies.
DEPS=(boost gsl netcdf4 openmm pyqt tbb tbb-devel)

# Where the Conda environment is stored.
CONDA_ENV=.conda_env

# Get the Sire version.
SIRE_VER=$(git --git-dir=$SRC_DIR/.git describe --tags | tr - _)

# Store the conda environment.
$HOME/sire.app/bin/conda env export -n base > $CONDA_ENV

# Loop over all dependences and replace with the version installed
# within the Conda enviroment.
echo "Updating Python dependencies..."
for dep in ${DEPS[@]}; do
    ver=$(grep "\- $dep=" $CONDA_ENV | awk -F "=" '{print $2}')
    sed -i.bak -e  "0,/$dep/s//$dep $ver/" -- $RECIPE && rm -- $RECIPE.bak
    echo "  $dep $ver"
done

# Update the Sire version number.
echo "Updating Sire version number: '$SIRE_VER'"
sed -i.bak "s/VERSION/$SIRE_VER/" -- $RECIPE && rm -- $RECIPE.bak

# Remove the Conda environment file.
rm -f .conda_env

# Update the sha256 checksums.
echo "Updating package checksums..."
# Check the OS.
if [ "$(uname)" == "Darwin" ]; then
    CHECKSUM_OSX=$(openssl sha256 $ARCHIVE_OSX | awk '{print $2}')
    echo "  macOS: $CHECKSUM_OSX"
    sed -i.bak -e "s/SHA256_OSX/$CHECKSUM_OSX/" -- $RECIPE && rm -- $RECIPE.bak
else
    CHECKSUM_LINUX=$(openssl sha256 $ARCHIVE_LINUX | awk '{print $2}')
    echo "  Linux: $CHECKSUM_LINUX"
    sed -i.bak -e "s/SHA256_LINUX/$CHECKSUM_LINUX/" -- $RECIPE && rm -- $RECIPE.bak
fi

echo "Recipe updated!"
