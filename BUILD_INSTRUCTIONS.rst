Sire build instructions
***********************

The following page details the Sire build process. Sire is built using
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__, creating
an isolated Conda environment into which Sire with be compiled and deployed.
This avoids the need for external dependencies, since everything that is
required is already available within the Conda ecosystem.

The main script used to compile Sire is
`compile_sire.sh <https://github.com/michellab/Sire/blob/devel/compile_sire.sh>`__
on Linux and `compile_sire.sh <https://github.com/michellab/Sire/blob/devel/compile_sire.bat>`__
on Windows. You shouldn't need to modify this file, other than updating
the ``MINICONDA_VERSION`` variable if needed. (Note that since we use
a specific Miniconda we only currently build Sire for Python version 3.7.)

The `build_sire.py <https://github.com/michellab/Sire/blob/devel/build/build_sire.py>`__
Python script is called by the above shell scripts to perform the build.
It is within this file that external dependencies are installed, using Conda
where possible, and pip if a package is not yet available via the Anaconda
Cloud. Please use the `conda-forge <https://conda-forge.org>`__ channel
as your first port of call, since this will help maximise compatibility between
packages. (We plan on eventually migrating to the conda-forge build platform
ourselves, so it is important that we can eventually depend on a single Conda
channel.) When choosing / updating package versions please choose the most
recent version that is supported on all operating systems that we build for,
i.e. ``linux-64``, ``osx-64``, and ``win-64``. To check which versions are
compatible use the Anaconda Cloud, e.g. for the `boost <https://anaconda.org/conda-forge/boost>`__.
(Obviously you can use an older version of a package if the latest is broken,
or if dependencies can't be resolved.) Note that versions aren't just
important for Sire itself, but also any applications that are built on top
of Sire, such as `BioSimSpace <https://github.com/michellab/biosimspace>`__.
Please don't change version numbers without checking that applications that
depend on them still work. (For example, `RDKit <https://www.rdkit.org>`__
is used by BioSimSpace, which depends on some of the same packages used
by Sire.)

We use `Azure Pipelines <https://dev.azure.com/michellab/Sire/_build>`__ for
continuous integration and deployment. On Linux this is manage using
Docker containers, with details for each stage of the build available
in the `docker <https://github.com/michellab/Sire/tree/devel/docker>`__
directory. Note that the pipeline can be unreliable because of network
timeouts, so if a build fails it is often worth trying to queue another
build of the same commit from the online control panel. (Check any error
messages first to see if it looks like a genuine issue.)

The `sire-conda <https://github.com/michellab/Sire/tree/devel/docker/sire-conda>`__
directory contains all of the files that are used to build and deploy
the Sire `Conda package <https://anaconda.org/michellab/sire>`__ during
the Azure Pipeline. The scripts will query the Sire Miniconda for the
versions of any external dependencies, then add them to the packages
*recipe*. This means that any version changes that you make in the
`build_sire.py <https://github.com/michellab/Sire/blob/devel/build/build_sire.py>`__
file will automatically be updated in the Conda package, i.e. you
only need to maintain dependency versions in a single location.
(The only thing that is needed is to update the list of external Conda
dependences in `update_recipe.sh <https://github.com/michellab/Sire/blob/devel/docker/sire-conda/update_recipe.sh>`__
if a new dependency is added.) Similarly, the Conda package version number is
also automatically inferred from the git commit (we use the tag as the version
and the number of commits since the tag as the build number). The Sire Conda
package is created by extracting all of the compiled Sire binary and libray
files from the Miniconda installation, then uploading them to the Oracle Cloud
as an archive. These are then downloaded and unpacked into the correct location
during the Conda installation. See `this <https://github.com/michellab/Sire/blob/devel/docker/sire-conda/create_package_file.sh>`__
script to see how the archive of files is created and `this <https://github.com/michellab/Sire/blob/devel/docker/sire-conda/recipe/build.sh>`__ script to see how the same files are unpacked in place
during the Conda install.
