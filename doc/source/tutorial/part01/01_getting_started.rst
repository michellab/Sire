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

   >>> mols = sr.load("https://siremol.org/m/cholesterol.sdf")
   >>> print(mols)
   XXXX

Molecules are loaded into a :class:`~Sire.System.System`. You can see how
many molecules have been loaded using the `.nMolecules()` function;

.. code-block:: python

   >>> print(mols.nMolecules())
   1

In this case, one molecule has been loaded. You can access this molecule via;

.. code-block:: python

   >>> mol = mols[0]
   >>> print(mol)
   XXXX

There are many ways to view the atoms in the molecule. One is to use
the index, e.g.

.. code-block:: python

   >>> atom = mol[0]
   >>> print(atom)
   XXXX

or

.. code-block:: python

   >>> atom = mol.atom(0)
   >>> print(atom)
   XXXX

would access the first atom in the molecule. The `.nAtoms()`
function returns the total number of atoms.

.. code-block:: python

   >>> print(mol.nAtoms())
   XXXX

You can loop over all of the atoms via the `.atoms()` function e.g.

.. code-block:: python

   >>> for atom in mol.atoms():
   ...     print(atom)
   XXXX

Molecules can be divided into residues, chains and segments. A residue
is a collection of atoms, a chain is a collection of residues, and a segment
is an arbitrary, but often-larger collection of atoms within a molecule.

You can access residues, chains and segments in similar ways to accessing
atoms, e.g.

.. code-block:: python

   >>> res = mol.residue(0)
   >>> print(res)
   XXXX

   >>> for res in mol.residues():
   ...     print(res)

You access atoms in a residue, chain or segment in a similar way, e.g.

.. code-block:: python

   >>> res = mol.residue(0)
   >>> atom = res.atom(0)
   >>> print(atom)
   XXXX

   >>> for atom in res.atoms():
   ...     print(atom)
   XXXX

Saving a molecule
-----------------

You save molecules using the :func:`Sire.save` function;

.. code-block:: python

   >>> sr.save(mol, "cholesterol.pdb")
   XXXX

Sire will automatically try to guess the file format from the file
extension. In this case, the molecule is saved in PDB format.

You can specify the format using the `format` argument.

.. code-block:: python

   >>> sr.save(mol, "cholesterol", format="mol2")
   XXXX

Note how the file format extension has been added automatically, and
that the full path to the file that was written is returned.

Loading from multiple files
---------------------------

It is often the case that molecular information needs to be read from
multiple files, e.g. a separate topology and coordinate file.

You load from multiple files simply by passing multiple filenames and/or
URLs to :func:`Sire.load`.

.. code-block:: python

   >>> mols = sr.load("https://siremol.org/m/urea.top",
   ...                "https://siremol.org/m/urea.gro")
   >>> print(mols)
   XXXX

You can pass in the filenames as multiple arguments or as a list,
whichever you find easiest.

.. code-block:: python

   >>> mols = sr.load(["https://siremol.org/m/urea.top",
   ...                 "https://siremol.org/m/urea.gro"])




Saving to multiple files
------------------------


Symmetric Input / Output
------------------------

One of Sire's design principles is that molecular file input and output
is symmetrical. This means that Sire can read in and write out the same
amount of information from a file (i.e. it can always read what it writes).

Another design principle is that information should not be lost. As much
as possible, Sire will load and preserve all information it can
read from a molecular file.

