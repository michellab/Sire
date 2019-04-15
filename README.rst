****
`Sire <http://siremol.org>`__
****

.. image:: https://dev.azure.com/michellab/Sire/_apis/build/status/michellab.Sire?branchName=devel
   :target: https://dev.azure.com/michellab/Sire/_build

.. image:: https://anaconda.org/michellab/sire/badges/downloads.svg
   :target: https://anaconda.org/michellab/sire

About
=====
`Sire <http://siremol.org>`__ is a free, open source, multiscale
molecular simulation framework, written to allow computational
modellers to quickly prototype and develop new algorithms for
molecular simulation and molecular design. Sire is written
as a collection of libraries, each of which contains self-contained
and robust C++/Python building blocks. These building blocks are
vectorised and thread-aware and can be streamed (saved/loaded)
to and from a version-controlled and tagged binary format,
thereby allowing them to be combined together easily to build
custom multi-processor molecular simulation applications.

For more information about how to use Sire, and about application
built with Sire, please `visit the Sire website <http://siremol.org>`__.

Installation
============

The easiest way to install Sire is using our `conda channel <https://anaconda.org/michellab/repo>`__::

    conda install -c conda-forge -c omnia -c michellab sire

To install the latest development version you can use::

    conda install -c conda-forge -c omnia -c michellab/label/dev sire

Unless you add the required channels to your Conda configuration, then you'll
need to add them when updating, e.g., for the development package::

    conda update -c conda-forge -c omnia -c michellab/label/dev sire

There are also many `pre-built binary packages <http://siremol.org/pages/binaries.html>`__,
which are available for Linux, Mac OS X and Windows, which are quick and easy to install.

However, as you are here, it is likely you want to download the latest,
greatest version of the code, which you will need to compile. To compile Sire,
you need a working C++ compiler with at least C++ 2014 support (gcc >= 5 or clang >= 3.7),
`cmake <http://cmake.org>`__
(version 3.0.0 or above), a Git client to download the source,
and a working internet connection (needed by
the Sire compilation scripts to download additional dependencies).

The easy install option is::

    git clone git@github.com:michellab/Sire.git
    cd Sire
    ./compile_sire.sh

A small word of warning, the compilation can easily take over an hour!

The above will download and install a new miniconda python installation,
into which Sire will be compiled and deployed (together with its
dependencies). This is by far the easiest way to compile and install Sire,
and is the route we strongly recommend. If you have any problems with
compiling and installing Sire, then please get in touch using the links below.

If you want to install Sire into an existing miniconda or Anaconda
Python installation, please follow the instructions in `build/INSTALL_INTO_ANACONDA.rst`.

Docker images
=============

If you don't want to build or install, you can also run Sire via one of our
docker images. The easy way to run the latest development image of Sire is via::

    docker run -it siremol/sire-devel:latest

This will download the latest Sire development container, and will run it,
giving you a bash prompt inside the container.

OpenMM compatibility
====================

Some Sire functionality requires `OpenMM <http://openmm.org>`__. Although
a bundled version is provided as part of the installation, this may not
be appropriate for your GPU drivers. To automatically detect and install
a suitable version of OpenMM, simply run the following command post-install::

    optimise_openmm

(Note that, depending on your installation method, `optimise_openmm` may
be location in `$HOME/sire.app/bin`.)

Support and Development
=======================

Bugs, Comments, Questions
--------------------------
For bug reports/sugguestions/complains please file an issue on
`GitHub <http://github.com/michellab/Sire>`__.
or contact the developers via the google user group: `https://groups.google.com/forum/#!forum/sire-users`

Developers guide
-----------------
Please `visit the website <http://siremol.org>`__ for information on how to
develop applications using Sire.

Azure Pipelines -- Autobuild feature
---------------------------
Since Sire is quite large, a build can take quite long and might not be neccessary
if a commit is only fixing a couple of typos. Simply add the line `***NO_CI***`
to your commit message and Azure Pipelines will not invoke an autobuild.

Note that every time you commit to devel, it will trigger a build of Sire,
full testing, construction of a package and upload to siremol.org (so that it
can be downloaded as the latest version of sire_devel_latest_linux.run). Please
think twice before committing directly to devel. You should ideally be working
in a feature branch, and only commit to devel once you are happy the code
works on your branch. Use `***NO_CI***` until you are happy that you want to
trigger a full build, test and deployment. This full pipeline will take
about 30 minutes to complete.

Have fun :-)
