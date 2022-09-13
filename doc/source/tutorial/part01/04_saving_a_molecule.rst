=================
Saving a molecule
=================

You save molecules using the :func:`sire.save` function;

>>> sr.save(mol, "cholesterol.pdb")
['/path/to/cholesterol.pdb']

``sire`` will automatically try to guess the file format from the file
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

Saving to strings
=================

You can also save a molecule to memory using the
:func:`sire.save_to_string` function. This saves the molecules
to memory, returning the lines of the file as a Python list.

You need to specify the format you want to use, e.g.
``pdb`` for PDB, ``gro`` for Gro87 etc.

>>> sr.save_to_string(mols, format="pdb")
['MODEL     1',
 'ATOM      1  C   MOL     1      -0.019   1.525   0.010  1.00  0.00           C',
 'ATOM      2  C   MOL     1       0.002  -0.004   0.002  1.00  0.00           C',
 'ATOM      3  C   MOL     1       0.658  -0.518   1.285  1.00  0.00           C',
 'ATOM      4  C   MOL     1       1.418   2.055  -0.000  1.00  0.00           C',
 'ATOM      5  C   MOL     1       2.164   1.442   1.173  1.00  0.00           C',
 ...
 'ATOM     72  H   MOL     1       3.900  -0.359   0.167  1.00  0.00           H',
 'ATOM     73  H   MOL     1       2.811  -1.766   0.235  1.00  0.00           H',
 'ATOM     74  H   MOL     1      -0.762   2.962  -1.208  1.00  0.00           H',
 'ENDMDL',
 'END']

>>> sr.save_to_string(mols, format="gro")
'cholesterol',
 '   74',
 '    1MOL      C    1  -0.001870   0.152540   0.001040',
 '    1MOL      C    2   0.000210  -0.000410   0.000200',
 '    1MOL      C    3   0.065850  -0.051780   0.128510',
 '    1MOL      C    4   0.141830   0.205520  -0.000040',
 '    1MOL      C    5   0.216410   0.144230   0.117300',
 '    1MOL      C    6   0.280940   0.223720   0.197530',
...
 '    1MOL      H   72   0.389970  -0.035920   0.016740',
 '    1MOL      H   73   0.281060  -0.176550   0.023460',
 '    1MOL      H   74  -0.076160   0.296250  -0.120770',
 '   0.00000   0.00000   0.00000']
