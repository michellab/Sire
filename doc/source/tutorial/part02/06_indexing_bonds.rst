==============
Indexing Bonds
==============

Bonds represent the chemical bonds between atoms in a molecule. A
:class:`~sire.mol.Bond` is a molecular container that contains the
two atoms that make up the bond.

For example, let's look at the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]
>>> print(mol)

