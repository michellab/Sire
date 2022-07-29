=============================
Part 3 - Molecular Properties
=============================

The :class:`~sire.mol.MoleculeView`-derived classes, such as
:class:`~sire.mol.Atom`, :class:`~sire.mol.Residue` etc., are
containers for molecular information. This information is held
via properties that are associated with each view in the molecule.

To see how this works, we will first load up the `aladip` system.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

Properties are accessed via the `property` function. This takes,
as argument, the name of the property you want to retrieve.
For example, the coordinates are held in the `coordinates` property.

>>> print(mol.property("coordinates"))
...

The coordinates

.. toctree::
   :maxdepth: 1

   part03/01_atom_properties
   part03/02_residue_properties
   part03/03_bond_properties
