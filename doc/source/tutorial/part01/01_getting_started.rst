===============
Getting Started
===============

Sire is Python library that is designed to make it easier for you
to build and manipulate molecular systems.

This tutorial assumes that you have installed ``Sire`` following
:doc:`the instructions here. <../../install>`, and have activated
the Anaconda / miniconda Python environment into which Sire was
installed.

.. warning::

  This tutorial is written for Sire version |SireVersion| or
  higher. If you are using an older version then please upgrade.

Importing Sire
--------------

Start a Python session (e.g. an interactive python console or a
Jupyter notebook).

You import Sire by typing

.. code-block:: python

   >>> import Sire

.. note::

   Note how this tutorial uses ``>>>`` to show a command that you should
   type into an interactive Python session, such as an ipython console or
   Jupyter notebook

As a convention, we will import the :mod:`Sire` Python module under the alias
`sr` to reduce the amount of typing. We will use `sr` throughout this tutorial to
mean :mod:`Sire`.

.. code-block:: python

   >>> import Sire as sr

Loading a molecule
------------------

We load molecules using the :func:`Sire.load` function. This accepts either
a filename, a URL, or a `PDB code <https://www.rcsb.org>`__.

For example, let's load a cholesterol molecule from
`https://siremol.org/m/cholesterol.sdf <https://siremol.org/m/cholesterol.sdf>`__.

.. code-block:: python

   mols = sr.load("https://siremol.org/m/cholesterol.sdf")

   print(mols)

This should load the molecule and output something like;

..

   Some output

Supports lots of file formats.

Saving a molecule
-----------------


