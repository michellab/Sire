#!/usr/bin/env bash

# Get the GitHub token and email.
GITHUB_TOKEN=${GITHUB_TOKEN:-$1}
GITHUB_EMAIL=${GITHUB_EMAIL:-$2}

# Set the Conda Forge feedstock directory.
CONDA_DIR=$HOME/staged-recipes

# Delete any existing Conda Forge directory.
if [ -d $CONDA_DIR ]; then
    echo "Deleting existing Conda directory: '$CONDA_DIR'"
    rm -rf $CONDA_DIR
fi

# Store the name of the recipe and template YAML files.
RECIPE=$CONDA_DIR/recipes/sire/meta.yaml
TEMPLATE=$CONDA_DIR/recipes/sire/template.yaml

# Clone the feedstock repository.
echo "Cloning Conda Forge feedstock: 'https://github.com/michellab/staged-recipes.git'"
git clone --branch devel https://github.com/michellab/staged-recipes.git $CONDA_DIR > /dev/null 2>&1

# Overwite the recipe with the template file.
cp $TEMPLATE $RECIPE

ARCHIVE_LINUX=$HOME/sire_conda_latest_linux.tar.bz2
ARCHIVE_OSX=$HOME/sire_conda_latest_osx.tar.bz2

# Download the macOS Conda package files. The macOS Azure build is much faster
# than the Linux one so these files should always be created by the time we
# get to this stage. If not, then the macOS Conda package will be one commit
# behind the Linux package.
echo "Downloading macOS Conda package file"
curl --silent --insecure --location https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/LSe3OL7yKxPq2d5BgBVrvWFWdFlMzBG4VKUbbMahXMU/n/chryswoods/b/sire_releases/o/sire_conda_latest_osx.tar.bz2 --output $HOME/sire_conda_latest_osx.tar.bz2

# List of python dependencies.
DEPS=(boost gsl netcdf4 openmm pyqt tbb tbb-dev)

# Where the Conda environment is stored.
CONDA_ENV=.conda_env

# Get the Sire version.
SIRE_VER=$(git --git-dir=$HOME/Sire/.git describe --tags)

# Store the conda environment.
$HOME/sire.app/bin/conda env export -n base > $CONDA_ENV

# Loop over all dependences and replace with the version installed
# within the Conda enviroment.
echo "Updating Python dependencies..."
for dep in ${DEPS[@]}; do
    ver=$(grep "\- $dep=" $CONDA_ENV | awk -F "=" '{print $2}')
    sed -i.bak -e "s/$dep/$dep $ver/" -- $RECIPE && rm -- $RECIPE.bak
    echo "  $dep $ver"
done

# Update the Sire version number.
echo "Updating Sire version number: '$SIRE_VER'"
sed -i.bak -e "s/VERSION/$SIRE_VER/" -- $RECIPE && rm -- $RECIPE.bak

# Remove the Conda environment file.
rm -f .conda_env

# Update the sha256 checksums.
CHECKSUM_LINUX=$(openssl sha256 $ARCHIVE_LINUX | awk '{print $2}')
CHECKSUM_OSX=$(openssl sha256 $ARCHIVE_OSX | awk '{print $2}')
echo "Updating package checksums..."
echo "  Linux: $CHECKSUM_LINUX"
echo "  macOS: $CHECKSUM_OSX"
sed -i.bak -e "s/SHA256_LINUX/$CHECKSUM_LINUX/" -- $RECIPE && rm -- $RECIPE.bak
sed -i.bak -e "s/SHA256_OSX/$CHECKSUM_OSX/" -- $RECIPE && rm -- $RECIPE.bak

# Remove the recipe backup file.
rm -f $RECIPE.bak

# Change to Conda package directory and update git config.
cd $CONDA_DIR
git config user.name "BioSimSpaceBot" > /dev/null 2>&1
git config user.email "$GITHUB_EMAIL" > /dev/null 2>&1

# Commit the changes to the Conda recipe.
echo "Commiting changes to feedstock"
git commit -a -m "Updating Conda recipe." > /dev/null 2>&1
echo "Pushing changes to remote"
git push --repo https://biosimspacebot:$GITHUB_TOKEN@github.com/michellab/staged-recipes.git origin devel > /dev/null 2>&1
echo "Feedstock updated!"
