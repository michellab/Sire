****
`Sire <http://siremol.org>`__
****

.. image:: https://github.com/michellab/Sire/workflows/Build/badge.svg
   :target: https://github.com/michellab/Sire/actions?query=workflow%3ABuild)
   :alt: Build status

.. image:: https://anaconda.org/michellab/sire/badges/downloads.svg
   :target: https://anaconda.org/michellab/sire
   :alt: Downloads

.. image:: https://img.shields.io/badge/License-GPL%20v2-blue.svg
   :target: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
   :alt: License

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

The easiest way to install Sire is using our `conda channel <https://anaconda.org/michellab/repo>`__.
Sire is built using dependencies from `conda-forge <https://conda-forge.org/>`__,
so please ensure that the channel takes strict priority. We recommend using
`Miniforge <https://github.com/conda-forge/miniforge>`__.

To create a new environment:

.. code-block:: bash

    conda create -n sire -c conda-forge -c michellab sire
    conda activate sire

To install the latest development version you can use:

.. code-block:: bash

    conda create -n sire-dev -c conda-forge -c michellab/label/dev sire
    conda activate sire-dev

If you find that Conda is particularly slow to install or upgrade,
then we advise using `mamba <https://github.com/TheSnakePit/mamba>`__:

.. code-block:: bash

    conda install -c conda-forge mamba

You can then replace all ``conda`` commands with ``mamba``, e.g.:

.. code-block:: bash

    mamba create -n sire -c conda-forge -c michellab sire

However, as you are here, it is likely you want to download the latest,
greatest version of the code, which you will need to compile. To compile Sire,
you need a Git client to download the source, and a working internet connection
(needed by the Sire compilation scripts to download additional dependencies).

First, you need to create and activate a conda environment, e.g.

.. code-block:: bash

    conda create -n sire-dev
    conda activate sire-dev

We find that Conda is particularly slow to install or upgrade,
so we advise installing `mamba <https://github.com/TheSnakePit/mamba>`__:

.. code-block:: bash

    conda install -c conda-forge mamba

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

    git clone git@github.com:michellab/Sire.git
    cd Sire
    python setup.py install

A small word of warning, the compilation can easily take over an hour!

The above will compile Sire in your existing conda environment.

If you plan to build `BioSimSpace <https://github.com/michellab/BioSimSpace>`__
on top of Sire, then you will need to resolve BioSimSpace's dependencies at
the time Sire is installed to ensure that it is built in a self-consistent way.
This can be achieved as follows:

.. code-block:: bash

    python setup.py --install-bss-deps install

Support and Development
=======================

Bugs, Comments, Questions
-------------------------
For bug reports/sugguestions/complains please file an issue on
`GitHub <http://github.com/michellab/Sire/issues>`__.
or contact the developers via the google user group: `https://groups.google.com/forum/#!forum/sire-users`

Developers guide
----------------
Please `visit the website <http://siremol.org>`__ for information on how to
develop applications using Sire.

GitHub actions
--------------
Since Sire is quite large, a build can take quite long and might not be neccessary
if a commit is only fixing a couple of typos. Simply add ``ci skip``
to your commit message and GitHub actions will not invoke an autobuild.

Note that every time you commit to devel, it will trigger a build of Sire,
full testing, construction of a Conda package and upload to our Anaconda
channel. Please think twice before committing directly to devel. You should
ideally be working in a _feature_ branch, and only commit to devel once you are
happy the code works on your branch. Use ``ci skip`` until you are happy that
you want to trigger a full build, test and deployment. This full pipeline will
take several hours to complete.

Have fun :-)
