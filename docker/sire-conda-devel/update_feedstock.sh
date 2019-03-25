#!/usr/bin/env bash

# Switch to the home directory.

# Get the GitHub token and email.
GITHUB_TOKEN=${GITHUB_TOKEN:-$1}
GITHUB_EMAIL=${GITHUB_EMAIL:-$2}

# Set the Conda Forge feedstock directory.
CONDA_DIR=$HOME/staged-recipes

# Delete any existing Conda Forge directory.
if [ -d $CONDA_DIR ]; then
    rm -rf $CONDA_DIR
fi

# Store the name of the recipe and template yaml files.
RECIPE=$CONDA_DIR/recipes/sire/meta.yaml
TEMPLATE=$CONDA_DIR/recipes/sire/template.yaml

# MacOS.
if [ "$(uname)" = "Darwin" ]; then
    ARCHIVE=$HOME/sire_conda_latest_osx.tar.bz2
    SHA256=SHA256_OSX

# Linux.
elif [ "$(uname)" = "Linux" ]; then
    ARCHIVE=$HOME/sire_conda_latest_linux.tar.bz2
    SHA256=SHA256_LINUX

    # List of python dependencies.
    DEPS=(boost gsl netcdf4 openmm pyqt tbb tbb-dev)

    # Where the Conda environment is stored.
    CONDA_ENV=.conda_env

    # Get the Sire version.
    SIRE_VER=$(git --git-dir=$HOME/Sire/.git describe --tags)

    # Store the conda environment.
    $HOME/sire.app/bin/conda env export -n base > $CONDA_ENV

    # Clone the feedstock repository.
    git clone --branch devel https://github.com/michellab/staged-recipes.git $CONDA_DIR > /dev/null 2>&1

    # Overwite the recipe with the template file.
    cp $TEMPLATE $RECIPE

    # Loop over all dependences and replace with the version installed
    # within the Conda enviroment.
    for dep in ${DEPS[@]}; do
        ver=$(grep "\- $dep=" $CONDA_ENV | awk -F "=" '{print $2}')
        sed -i.bak -e "s/$dep/$dep $ver/" -- $RECIPE && rm -- $RECIPE.bak
    done

    # Update the Sire version number.
    sed -i.bak -e "s/VERSION/$SIRE_VER/" -- $RECIPE && rm -- $RECIPE.bak

    # Remove the Conda environment file.
    rm -f .conda_env

# Unsupported OS.
else
    echo "Unsupported operating system: '$uname'"
    exit -1
fi

# Update the sha256 checksum.
CHECKSUM=$(openssl sha256 $ARCHIVE | awk '{print $2}')
sed -i.bak -e "s/$SHA256/$CHECKSUM/" -- $RECIPE && rm -- $RECIPE.bak

# Remove the recipe backup file.
rm -f $RECIPE.bak

# Change to Conda package directory and update git config.
cd $CONDA_DIR
git config user.name "BioSimSpaceBot"
git config user.email "$GITHUB_EMAIL"

# Commit the changes to the Conda recipe.
git commit -a -m "Updating Conda recipe."
git push --repo https://biosimspacebot:$GITHUB_TOKEN@github.com/michellab/staged-recipes.git > /dev/null 2>&1
