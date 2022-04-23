=================
Saving a molecule
=================

You save molecules using the :func:`Sire.save` function;

>>> sr.save(mol, "cholesterol.pdb")
['/path/to/cholesterol.pdb']

Sire will automatically try to guess the file format from the file
extension. In this case, the molecule is saved in PDB format.

You can specify the format using the ``format`` argument.

>>> sr.save(mol, "cholesterol", format="mol2")
['/path/to/cholesterol.mol2']

Note how the file format extension has been added automatically, and
that the full path to the file that was written is returned.

You can specify multiple file formats, and thus write to multiple
files, e.g.

>>> sr.save(mol, "cholesterol", format=["mol2", "pdb"])
['/path/to/cholesterol.mol2', '/path/to/cholesterol.pdb']

or you can specify the filenames directly, e.g.

>>> sr.save(mol, ["chol.pdb", "chol.mol2"])
['/path/to/chol.pdb', '/path/to/chol.mol2']
