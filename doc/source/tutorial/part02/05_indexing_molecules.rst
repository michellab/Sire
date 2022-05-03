==================
Indexing Molecules
==================

A :class:`Sire.Mol.Molecule` holds all of the molecular information of
a single molecule. An input file may contain more than one molecule.
For example, the :class:`Sire.System.System` of molecules loaded from
these files...

>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
Downloading from 'https://siremol.org/m/ala.top'...
Downloading from 'https://siremol.org/m/ala.crd'...
>>> print(mols)
System( name=ACE num_molecules=631 num_residues=633 num_atoms=1912 )

contains 631 molecules.


