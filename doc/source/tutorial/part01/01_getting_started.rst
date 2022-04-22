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

Simple indexing
---------------

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

or

.. code-block:: python

   >>> for atom in mol.atoms(sr.range(0, 5)):
   ...     print(atom)

   Atom( C : 1 )
   Atom( C : 2 )
   Atom( C : 3 )
   Atom( C : 4 )
   Atom( C : 5 )

or

.. code-block:: python

   >>> for atom in mol.atoms([0, 2, 5, 8]):
   ...     print(atom)

   Atom( C : 1 )
   Atom( C : 3 )
   Atom( C : 6 )
   Atom( C : 9 )

.. note::

   The `sr.range` function is the standard Python range,
   except it returns a list rather than an iterator. This
   is needed because `mol.atoms()` can accept a list, but
   not yet a python iterator

Molecules can be divided into residues, chains and segments. A residue
is a collection of atoms, a chain is a collection of residues, and a segment
is an arbitrary, but often-larger collection of atoms within a molecule.

You can access residues, chains and segments in similar ways to accessing
atoms, e.g.

.. code-block:: python

   >>> res = mol.residue(0)
   >>> print(res)

   Residue( MOL : 1 )

   >>> for res in mol.residues():
   ...     print(res)

   Residue( MOL : 1 )

You access atoms in a residue, chain or segment in a similar way, e.g.

.. code-block:: python

   >>> res = mol.residue(0)
   >>> atom = res.atom(0)
   >>> print(atom)

   Atom( C : 1 )

   >>> for atom in res.atoms([0, 2, 4]):
   ...     print(atom)

   Atom( C : 1 )
   Atom( C : 3 )
   Atom( C : 5 )

Saving a molecule
-----------------

You save molecules using the :func:`Sire.save` function;

.. code-block:: python

   >>> sr.save(mol, "cholesterol.pdb")

   ['/path/to/cholesterol.pdb']

Sire will automatically try to guess the file format from the file
extension. In this case, the molecule is saved in PDB format.

You can specify the format using the `format` argument.

.. code-block:: python

   >>> sr.save(mol, "cholesterol", format="mol2")

   ['/path/to/cholesterol.mol2']

Note how the file format extension has been added automatically, and
that the full path to the file that was written is returned.

You can specify multiple file formats, and thus write to multiple
files, e.g.

.. code-block:: python

   >>> sr.save(mol, "cholesterol", format=["mol2", "pdb"])

   ['/path/to/cholesterol.mol2', '/path/to/cholesterol.pdb']

or you can specify the filenames directly, e.g.

.. code-block:: python

   >>> sr.save(mol, ["chol.pdb", "chol.mol2"])

   ['/path/to/chol.pdb', '/path/to/chol.mol2']

Loading from multiple files
---------------------------

It is often the case that molecular information needs to be read from
multiple files, e.g. a separate topology and coordinate file.

You load from multiple files simply by passing multiple filenames and/or
URLs to :func:`Sire.load`.

.. code-block:: python

   >>> mols = sr.load("https://siremol.org/m/ala.top",
   ...                "https://siremol.org/m/ala.crd")
   Downloading from 'https://siremol.org/m/ala.top'...
   Downloading from 'https://siremol.org/m/ala.crd'...

   >>> print(mols)

   System( name=ACE nMolecules=631 nResidues=633 nAtoms=1912 )

You can pass in the filenames as multiple arguments or as a list,
whichever you find easiest.

.. code-block:: python

   >>> mols = sr.load(["https://siremol.org/m/ala.top",
   ...                 "https://siremol.org/m/ala.crd"])


If the files or URLs have a common base, then you can save some typing
by using :func:`Sire.expand`, e.g.

.. code-block:: python

   >>> mols = sr.load(sr.expand("https://siremol.org/m",
   ...                          "ala.top", "ala.crd"))

or

.. code-block:: python

   >>> mols = sr.load(sr.expand(sr.tutorial_url,
   ...                          ["ala.top", "ala.crd"]))

.. note::

   `sr.tutorial_url` expands to the base URL for tutorial files
   (https://siremol.org/m). It is worth using this variable for
   the tutorial as it auto-completes and will reduce errors.

If you are loading files, you can also make use of glob expressions
(wildcard expansions), e.g.

.. code-block:: python

   >>> mols = sr.load("ala.*")

.. note::

   This line loads the `ala.top` and `ala.crd` files that
   were downloaded by the above lines. This is because Sire downloads
   files at URLs to the current directory. You can tell it to use
   a different directory by passing that in via the `directory`
   argument, e.g. `sr.load(sr.expand(sr.tutorial_url,"cholesterol.sdf"), directory="tmp")`.
   The directory will be created automatically if it doesn't exist.

Supported file formats
----------------------

Sire supports reading and writing to many common molecular file formats.
You can print the list of supported formats using

.. code-block:: python

   >>> print(sr.supported_formats())

   ## Parser Gro87 ##
   Supports files: gro
   Gromacs Gro87 structure format files.
   ##################

   ## Parser GroTop ##
   Supports files: top
   Gromacs Topology format files.
   ###################

   ## Parser MOL2 ##
   Supports files: mol2
   Sybyl Mol2 format files.
   #################

   ## Parser PDB ##
   Supports files: pdb
   Protein Data Bank (PDB) format files.
   ################

   ## Parser PRM7 ##
   Supports files: prm7, prm, top7, top, prmtop7, prmtop
   Amber topology/parameter format files supported from Amber 7 upwards.
   #################

   ## Parser PSF ##
   Supports files: psf
   Charmm PSF format files.
   ################

   ## Parser RST ##
   Supports files: rst, crd, trj, traj
   Amber coordinate/velocity binary (netcdf) restart/trajectory files supported since Amber 9, now default since Amber 16.
   ################

   ## Parser RST7 ##
   Supports files: rst7, rst, crd7, crd
   Amber coordinate/velocity text (ascii) restart files supported from Amber 7 upwards.
   #################

   ## Parser SDF ##
   Supports files: sdf, mol
   Structure Data File (SDF) format files.
   ################

   ## Parser SUPPLEMENTARY ##
   Supports files: *
   Files that are supplementary to a lead parser.
   ##########################

SHOULD ADD IN A COMMAND TO AUTOMATICALLY DOWNLOAD THE GROMACS FILES
IF GROMACS_HOME IS NOT SET

Symmetric Input / Output
------------------------

One of Sire's design principles is that molecular file input and output
is symmetrical. This means that Sire can read in and write out the same
amount of information from a file (i.e. it can always read what it writes).

Another design principle is that information should not be lost. As much
as possible, Sire will load and preserve all information it can
read from a molecular file.
