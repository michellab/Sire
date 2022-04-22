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

   Downloading from 'https://siremol.org/m/cholesterol.sdf'...

   >>> print(mols)

   System( name=cholesterol nMolecules=1 nResidues=1 nAtoms=74 )

Molecules are loaded into a :class:`~Sire.System.System`. You can see how
many molecules have been loaded using the `.nMolecules()` function;

.. code-block:: python

   >>> print(mols.nMolecules())

   1

In this case, one molecule has been loaded. You can access this molecule via;

.. code-block:: python

   >>> mol = mols[0]
   >>> print(mol)

   Molecule( 2.11 : nAtoms=74, nResidues=1 )

.. note::

   The `2.11` is a number that Sire uses to identify this molecule.
   We will explain what this number is and how it is formed in a
   later chapter. Note that your molecule may have a different
   identifier.

There are many ways to view the atoms in the molecule. One is to use
the index, e.g.

.. code-block:: python

   >>> atom = mol[0]
   >>> print(atom)

   Atom( C : 1 )

or

.. code-block:: python

   >>> atom = mol.atom(0)
   >>> print(atom)

   Atom( C : 1 )

would access the first atom in the molecule. The `.nAtoms()`
function returns the total number of atoms.

.. code-block:: python

   >>> print(mol.nAtoms())

   74

You can loop over all of the atoms via the `.atoms()` function e.g.

.. code-block:: python

   >>> for atom in mol.atoms():
   ...     print(atom)

   Atom( C : 1 )
   Atom( C : 2 )
   Atom( C : 3 )
   Atom( C : 4 )
   Atom( C : 5 )
   Atom( C : 6 )
   Atom( C : 7 )
   Atom( C : 8 )
   Atom( C : 9 )
   Atom( C : 10 )
   etc.

You can also loop over a slice of atoms, e.g.

.. code-block:: python

   >>> for atom in mol[0:5]:
   ...     print(atom)

   Atom( C : 1 )
   Atom( C : 2 )
   Atom( C : 3 )
   Atom( C : 4 )
   Atom( C : 5 )



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


If the files or URLs have a common base, then you can save some typing
by using :func:`Sire.expand`, e.g.

.. code-block:: python

   >>> mols = sr.load(sr.expand("https://siremol.org/m",
   ...                          ["urea.top", "urea.gro"]))

If you are loading files, you can also make use of glob expressions
(wildcard expansions), e.g.

.. code-block:: python

   >>> mols = sr.load("urea.*")

.. note::

   This line loads the `urea.top` and `urea.gro` files that
   were downloaded by the above lines. This is because Sire downloads
   files at URLs to the current directory. You can tell it to use
   a different directory by passing that in via the `directory`
   argument, e.g. `sr.load("cholesterol.sdf", directory="tmp")`.
   The directory will be created automatically if it doesn't exist.


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

