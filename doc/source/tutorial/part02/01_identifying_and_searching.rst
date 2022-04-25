======================
Indexing and Searching
======================

The core of Sire are the various :class:`~Sire.Mol.MoleculeView`-derived
classes, such as :class:`~Sire.Mol.Atom`, :class:`~Sire.Mol.Residue`,
:class:`~Sire.Mol.Chain`, :class:`~Sire.Mol.Segment` and
:class:`~Sire.Mol.Molecule`, amongst others.

These can all be considered as containers for molecular information.
:class:`~Sire.Mol.Atom` is a container for atomic information,
:class:`~Sire.Mol.Molecule` is a container for molecular information etc.

We access this information for indexing or searching into these containers.

For example, let's load up the protein ``7SA1`` from the `PDB <https://www.rcsb.org/structure/7SA1>`__

>>> import Sire as sr
>>> mols = sr.load("7SA1")
Downloading from 'https://files.rcsb.org/download/7SA1.pdb.gz'...
7SA1.pdb.gz
Unzipping './7SA1.pdb.gz'...
>>> mol = mols[0]
>>> print(mol)
Molecule( 2.10 : nAtoms=11728, nResidues=1518 )

.. note::

Sire automatically downloads and unpacks structures from the PDB. Just
put in the PDB code as the argument to :func:`Sire.load`.

Molecules are constructed as atoms, which be can be (optionally) arranged
into residues, chains and segments. We can get the number of each of
these using

>>> print(f"The number of atoms is {mol.nAtoms()}")
The number of atoms is 11728
>>> print(f"The number of residues is {mol.nResidues()}")
The number of residues is 1518
>>> print(f"The number of chains is {mol.nChains()}")
The number of chains is 4
>>> print(f"The number of segments is {mol.nSegments()}")
The number of segments is 0

.. note::

Unlike most Python functions, which are named using underscores,
Sire functions are named using camelCase.
So we have ``mol.nAtoms()``, not ``mol.num_atoms()``. We have
kept to using camelCase for C++ functions, and underscore_naming
for Python functions, so that you can know which are C++
(and thus fast and parallelisable), and which are Python.
