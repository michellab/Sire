****
`Sire <http://sire.openbiosim.org>`__
****

.. image:: https://github.com/OpenBioSim/sire/workflows/Build/badge.svg
   :target: https://github.com/OpenBioSim/sire/actions?query=workflow%3ABuild)
   :alt: Build status

.. image:: https://anaconda.org/OpenBioSim/sire/badges/downloads.svg
   :target: https://anaconda.org/OpenBioSim/sire
   :alt: Downloads

.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/old-licenses/gpl-3.0.en.html
   :alt: License

About
=====

Sire is a molecular modelling framework that provides extensive
functionality to manipulate representations of biomolecular systems.

It is used as a key component of `BioSimSpace <https://biosimspace.org>`__,
and is distributed and supported as an open source community project by
`OpenBioSim <https://openbiosim.org>`__.

For more information about how to use Sire, and about application
built with Sire, please `visit the Sire website <http://sire.openbiosim.org>`__.

Installation
============

The easiest way to install Sire is using our `conda channel <https://anaconda.org/openbiosim/repo>`__.
Sire is built using dependencies from `conda-forge <https://conda-forge.org/>`__,
so please ensure that the channel takes strict priority. We recommend using
`mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`__.

To create a new environment:

.. code-block:: bash

    mamba create -n openbiosim python==3.9
    mamba activate openbiosim
    mamba install -c openbiosim sire

To install the latest development version you can use:

.. code-block:: bash

    mamba create -n openbiosim-dev python==3.9
    mamba activate openbiosim-dev
    mamba install -c openbiosim/label/dev sire

However, as you are here, it is likely you want to download the latest,
greatest version of the code, which you will need to compile. To compile
sire,
you need a git client to download the source, and a working internet connection
(needed by the sire compilation scripts to download additional dependencies).

First, you need to create and activate a conda environment, e.g.

.. code-block:: bash

    mamba create -n openbiosim-dev python==3.9
    mamba activate openbiosim-dev

Next, you need to install the Sire build dependencies.

.. code-block:: bash

    mamba install cmake pip-requirements-parser

You will also need to install compilers, e.g. on Linux use

.. code-block:: bash

    mamba install gcc gxx

on MacOS use

.. code-block:: bash

    mamba install clang clangxx

and on Windows use

.. code-block:: bash

    mamba install conda-build

Next, you can clone the Sire source code and compile and install Sire::

    git clone https://github.com/OpenBioSim/sire
    cd sire
    python setup.py install

A small word of warning, the compilation can easily take over an hour!

The above will compile sire in your existing conda environment.

If you plan to build `BioSimSpace <https://github.com/michellab/BioSimSpace>`__
on top of sire, then you will need to resolve BioSimSpace's dependencies at
the time sire is installed to ensure that it is built in a self-consistent way.
This can be achieved as follows:

.. code-block:: bash

    python setup.py --install-bss-deps install

Support and Development
=======================

Bugs, Comments, Questions
-------------------------
For bug reports/sugguestions/complains please file an issue on
`GitHub <http://github.com/OpenBioSim/sire/issues>`__.

Developers guide
----------------
Please `visit the website <http://sire.openbiosim.org>`__ for information on how to
develop applications using sire.

GitHub actions
--------------
Since sire is quite large, a build can take quite long and might not be neccessary
if a commit is only fixing a couple of typos. Simply add ``ci skip``
to your commit message and GitHub actions will not invoke an autobuild.

Note that every time you commit to devel, it will trigger a build of sire,
full testing, construction of a Conda package and upload to our Anaconda
channel. Please think twice before committing directly to devel. You should
ideally be working in a _feature_ branch, and only commit to devel once you are
happy the code works on your branch. Use ``ci skip`` until you are happy that
you want to trigger a full build, test and deployment. This full pipeline will
take several hours to complete.

Have fun :-)
