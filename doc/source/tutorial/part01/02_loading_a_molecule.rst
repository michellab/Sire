==================
Loading a molecule
==================

We load molecules using the :func:`sire.load` function. This accepts either
a filename, a URL, or a `PDB code <https://www.rcsb.org>`__.

For example, let's load a cholesterol molecule from
`https://siremol.org/m/cholesterol.sdf <https://siremol.org/m/cholesterol.sdf>`__.

>>> mols = sr.load("https://siremol.org/m/cholesterol.sdf")
Downloading from 'https://siremol.org/m/cholesterol.sdf'...

>>> print(mols)
System( name=cholesterol num_molecules=1 num_residues=1 num_atoms=74 )

Molecules are loaded into a :class:`~sire.system.System`. You can see how
many molecules have been loaded using the :func:`~sire.mol.Molecule.num_molecules`
function;

>>> print(mols.num_molecules())
1

In this case, one molecule has been loaded. You can access this molecule via;

>>> mol = mols[0]
>>> print(mol)
Molecule( MOL:2   num_atoms=74 num_residues=1 )

This shows that the molecule is called ``MOL`` and has molecule number 2.
It contains 74 atoms in 1 residue.
