==================
Loading a molecule
==================

We load molecules using the :func:`Sire.load` function. This accepts either
a filename, a URL, or a `PDB code <https://www.rcsb.org>`__.

For example, let's load a cholesterol molecule from
`https://siremol.org/m/cholesterol.sdf <https://siremol.org/m/cholesterol.sdf>`__.

>>> mols = sr.load("https://siremol.org/m/cholesterol.sdf")
Downloading from 'https://siremol.org/m/cholesterol.sdf'...

>>> print(mols)
System( name=cholesterol nMolecules=1 nResidues=1 nAtoms=74 )

Molecules are loaded into a :class:`~Sire.System.System`. You can see how
many molecules have been loaded using the :func:`~Sire.Mol.Molecule.nMolecules`
function;

>>> print(mols.nMolecules())
1

In this case, one molecule has been loaded. You can access this molecule via;

>>> mol = mols[0]
>>> print(mol)
Molecule( 2.11 : nAtoms=74, nResidues=1 )

.. note::

   The ``2.11`` is a number that Sire uses to identify this molecule.
   We will explain what this number is and how it is formed in a
   later chapter. Note that your molecule may have a different
   identifier.
