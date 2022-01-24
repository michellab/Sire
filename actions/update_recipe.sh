#!/usr/bin/env bash

# Get the path to GitHub workspace.
SRC_DIR=${GITHUB_WORKSPACE:-$1}

# Set the Conda build directory path.
CONDA_DIR="$GITHUB_WORKSPACE"/recipes/sire

# Store the name of the recipe and template YAML files.
RECIPE="$CONDA_DIR"/meta.yaml
TEMPLATE="$CONDA_DIR"/template.yaml

# Overwite the recipe with the template file.
cp "$TEMPLATE" "$RECIPE"

# Get the Sire version. (Latest tag.)
SIRE_VERSION=$(git --git-dir="$SRC_DIR"/.git --work-tree="$SRC_DIR" describe --tags --abbrev=0)

# Get the build number. (Number of commits since last tag.)
SIRE_BUILD=$(git --git-dir="$SRC_DIR"/.git --work-tree="$SRC_DIR" log --oneline "$SIRE_VERSION".. | wc -l)

# Get the Sire branch.
SIRE_BRANCH=$(git rev-parse --abbrev-ref HEAD)

# Update the Sire version number.
echo "Updating Sire version number: '$SIRE_VERSION'"
sed -i.bak -e "s/VERSION/$SIRE_VERSION/" "$RECIPE" && rm "$RECIPE".bak

# Update the build number.
echo "Updating Sire build number: '$SIRE_BUILD'"
sed -i.bak -e "s/SIRE_BUILD/$SIRE_BUILD/" "$RECIPE" && rm "$RECIPE".bak

# Update the branch name.
echo "Updating Sire branch name: '$SIRE_BRANCH'"
sed -i.bak -e "s/SIRE_BRANCH/$SIRE_BRANCH/" "$RECIPE" && rm "$RECIPE".bak

echo "Recipe updated!"
