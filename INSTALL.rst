=========================
INSTALLATION INSTRUCTIONS
=========================

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

Packages are available for Python > 3.7 on Windows, MacOS and Linux.
