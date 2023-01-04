============
Installation
============

Binary packages for sire are available on MacOS, Linux and Windows.
Sire can be compiled on any UNIX or Windows-compatible operating system
running on X86-64, ARM64 or PowerPC processors.

You have a range of options for installing the software.

No-installation - Run in a Web Browser
======================================

We run a completely free `JupyterHub <https://try.openbiosim.org>`__ on
which we have sire installed.

This is at `try.openbiosim.org <https://try.openbiosim.org>`__.
You only need a `GitHub account <https://github.com>`__, which is
used to log into the server.

Simply go to `try.openbiosim.org <https://try.openbiosim.org>`__ in your
web browser, and log in using your `GitHub account <https://github.com>`__.
This will start a Jupyter Lab instance. In here you can start a terminal,
and then use ``sire`` directly via ``ipython``. Or you can start a Jupyter
notebook and use ``sire`` there.

To import ``sire``, at a Python prompt type

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use ``sire``.

.. note::

   This free JupyterHub server is limited. You only have up to 2 GB of
   memory and at most 1 processor core. Disk storage is temporary,
   and any data will be lost when you log out. Because it only
   supports a limited number of concurrent users, inactive sessions will be
   automatically stopped and logged out after 20 minutes. Please only
   use this service to explore and learn ``sire``.
   Do not use it for production work.

Easy installation - Run in a conda environment
==============================================

The easiest way to install ``sire`` is in a new
`conda environment <https://anaconda.org>`__.

You can use any conda environment or installation. We recommend using
`mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`__,
as this is pre-configured to use `conda-forge <https://conda-forge.org>`__,
and bundles `mamba <https://mamba.readthedocs.io/en/latest/>`__, which
is a fast drop-in replacement for `conda <https://conda.io>`__.

Option 1. Installing a new copy of ``mambaforge``
-------------------------------------------------

To install a new copy of
`mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`__,
first download a ``Mambaforge`` from
`this page <https://github.com/conda-forge/miniforge#mambaforge>`__ that
matches your operating system and processor.

Install ``Mambaforge`` following the
`instructions here <https://github.com/conda-forge/miniforge#install>`__.

Once installed, you should be able to run the ``mamba`` command to
install other packages (e.g. ``mamba -h`` will print out help on
how to use the ``mamba`` command).

Option 2. Using an existing anaconda/miniconda install
------------------------------------------------------

If you want to use an existing anaconda or miniconda installation,
then first open a terminal with that distribution activated.
For example, open a terminal via anaconda navigator, or
open a terminal and run
``source /path/to/conda/bin/activate``, where ``/path/to/conda`` is
the full path to your anaconda or miniconda installation.

You should now be able to run the ``conda`` command to install other
packages (e.g. ``conda -h`` will print out help on how to use the
``conda`` command). We highly recommend that you use ``mamba`` as a
drop-in replacement for ``conda``, so first install ``mamba``.

.. code-block:: bash

   $ conda install -c conda-forge mamba

This should install mamba. If this fails, then your anaconda or miniconda
environment is likely quite full, or else it is outdated. We recommend
going back to Option 1 and installing a new copy of ``mambaforge``.

If this works, then you should now be able to run the ``mamba`` command
to install other packages (e.g. ``mamba -h`` will print out help
on how to use the ``mamba`` command).

Installing sire into a new environment
--------------------------------------

We recommend that ``sire`` is installed into a new (clean) environment.
This minimises the risk of failures caused by incompatible dependencies.

Sire is currently packaged for Python 3.8 and Python 3.9. We will start
by creating a Python 3.9 environment that we will call ``openbiosim``.

.. code-block:: bash

   $ mamba create -n openbiosim python==3.9

We can now install ``sire`` into that environment by typing

.. code-block:: bash

   $ mamba install -n openbiosim -c openbiosim sire

.. note::

   The option ``-n openbiosim`` tells ``mamba`` to install ``sire``
   into the ``openbiosim`` environment. The option ``-c openbiosim``
   tells ``mamba`` to install ``sire`` from the ``openbiosim``
   conda channel.

You may (optionally) want to install additional tools such as
``ipython`` and ``jupyterlab``. To do this, type

.. code-block:: bash

   $ mamba install -n openbiosim ipython jupyterlab

To run ``sire``, you must now activate the ``openbiosim`` environment.
You can do this by typing

.. code-block:: bash

   $ conda activate openbiosim

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import ``sire`` by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use ``sire``.

Also easy installation - Run in a container
===========================================

Another route to install ``sire`` is to download and run our
pre-built containers. These can be run via
`docker <https://www.docker.com>`__ (on Linux, MacOS and Windows)
or via `podman <https://podman.io>`__ (on Linux).

To run via `docker <https://www.docker.com>`__, simply type;

.. code-block:: bash

   $ docker run -p 8888:8888 -it openbiosim/sire:latest

or, via `podman <https://podman.io>`__, type;

.. code-block:: bash

   $ podman run -p 8888:8888 -it openbiosim/sire:latest

This will download the container from
`hub.docker.com <https://hub.docker.com/r/openbiosim/sire>`__ and
will start a command prompt in that container.

You can now type ``python``, ``ipython`` or ``jupyter lab``
to start a python, ipython or jupyter lab session.

.. note::

   The option ``-p 8888:8888`` tells docker/podman to redirect
   port ``8888`` on your computer to port ``8888`` in the
   container. This will let you open a browser and navigate to
   the URL printed by ``jupyter lab`` if you are using jupyter.
   You can drop this option if you don't want to use
   ``jupyter lab``.

.. note::

   You can map directories from your computer into the container
   by using the ``-v`` option. For example,
   ``-v $HOME/input:/home/openbiosim/input`` would map your
   ``input`` folder in your home directory to the ``input`` folder
   in the home directory of the container. This will let ``sire``
   read and write files on your computer.

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import ``sire`` by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use ``sire``.

Harder installation - Compile from source
=========================================

Sometimes you will want to compile and run ``sire`` from source.
This could be because we don't distribute a binary package for your
operating system, or because you want to use a newer version
(e.g. code from the ``devel`` branch, or from your own feature
branch if you are a developer).

You compile ``sire`` into an existing anaconda / miniconda environment.
Please create and activate an environment, e.g. by following
the "Option 1" instructions to install a fresh ``mambaforge`` and
then creating and activating Python 3.9 environment called
``openbiosim``.

Next, download the source code. You could download the latest development
version of ``sire`` by typing;

.. code-block:: bash

   $ git clone https://github.com/openbiosim/sire

This will download into a directory called ``sire``. Navigate into
this directory (e.g. ``cd sire``).

.. note::

   This will fail if ``git`` is not installed on your computer.
   You can easily install ``git`` using ``mamba``, e.g.
   run ``mamba install git``.

You can change to a different branch using the ``git checkout BRANCH``
command, e.g.

.. code-block:: bash

   $ git checkout main

will check out the ``main`` branch of ``sire``. This always corresponds
to the last released version of ``sire``. Or, you can check out a
feature branch using

.. code-block:: bash

   $ git checkout feat_name

where ``feat_name`` should be replaced by the name of the feature
branch you want to compile.

Compilation and installation of ``sire`` is managed via the
`setup.py <https://github.com/michellab/Sire/blob/devel/setup.py>`__
script.

Run

.. code-block:: bash

   $ python setup.py --help

to get a help on all of the options.

Typically, you just want to compile and install ``sire``. To do this,
type

.. code-block:: bash

   $ python setup.py install

This will download and install all of the dependencies via ``mamba``
(or ``conda`` if you haven't installed ``mamba``). It will then compile
the ``sire`` C++ libraries, and then the Python wrappers. Be patient,
as compilation can take quite a while!

.. note::

   You need to have Visual Studio 2017 C++ compiler installed to compile on Windows.
   The easiest way to do this is to `install chocolatey <https://chocolatey.org/install>`__
   and then install the compilers using the command 
   ``choco install visualstudio2017-workload-vctools``. This is all free, but 
   you will need admin access to install chocolatey.

If you plan to install `BioSimSpace <https://biosimspace.org>`__ on
top of ``sire``, then you should install using;

.. code-block:: bash

   $ python --install-bss-deps install

This will use ``mamba`` (or ``conda``) to download and install all of
BioSimSpace's dependencies as well. This ensures that incompatible versions
of shared dependencies are not accidentally installed.

Once ``sire`` has installed, you can import it in a ``python``,
``ipython`` or ``jupyter lab`` session by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use ``sire``.
