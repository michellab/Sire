=================
Indexing Residues
=================

Residues are collections of atoms. They typically represent an amino
acid residue in a protein. Residues are implemented via the
:class:`~Sire.Mol.Residue` class, which itself is a molecular container
for :class:`~Sire.Mol.Atom` objects.

You can access residues in a molecule container using the
:func:`~Sire.Mol.Residue.residue` and :func:`~Sire.Mol.Residue.residues`
functions.

>>> print(mol.residue(0))
Residue( ILE:6   nAtoms=8 )

gives the molecule at index 0, while

>>> print(mol.residue("ALA"))
Selector<SireMol::Residue>( size=155
0:  Residue( ALA:23  nAtoms=5 )
1:  Residue( ALA:30  nAtoms=5 )
2:  Residue( ALA:53  nAtoms=5 )
3:  Residue( ALA:65  nAtoms=5 )
4:  Residue( ALA:85  nAtoms=5 )
...
150:  Residue( ALA:578 nAtoms=5 )
151:  Residue( ALA:584 nAtoms=5 )
152:  Residue( ALA:593 nAtoms=5 )
153:  Residue( ALA:646 nAtoms=5 )
154:  Residue( ALA:691 nAtoms=5 )
)

returns all residues that are named "ALA".

The :func:`~Sire.Mol.Residue.residue` will raise a KeyError if more than
one residue matches the search.

>>> print(mol.residue("ALA"))
---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
Input In [33], in <cell line: 1>()
----> 1 print(mol.residue("ALA"))
<BLANKLINE>
KeyError: "SireMol::duplicate_residue: More than one residue matches the ID ResName('ALA') (number of matches is 155). (call Sire.Error.get_last_error_details() for more info)"

You can slice residues using the ``range`` function, e.g.

>>> print(mol.residues(range(0, 10)))
Selector<SireMol::Residue>( size=10
0:  Residue( ILE:6   nAtoms=8 )
1:  Residue( VAL:7   nAtoms=7 )
2:  Residue( LEU:8   nAtoms=8 )
3:  Residue( LYS:9   nAtoms=9 )
4:  Residue( SER:10  nAtoms=6 )
5:  Residue( SER:11  nAtoms=6 )
6:  Residue( ASP:12  nAtoms=8 )
7:  Residue( GLY:13  nAtoms=4 )
8:  Residue( VAL:22  nAtoms=7 )
9:  Residue( ALA:23  nAtoms=5 )
)

The result, a :class:`~Sire.Mol.Selector_Residue_` is also a molecular
container, and can be used like :class:`~Sire.Mol.Selector_Atom_`.

>>> print(mol.residues("ALA")[0:5])
Selector<SireMol::Residue>( size=5
0:  Residue( ALA:23  nAtoms=5 )
1:  Residue( ALA:30  nAtoms=5 )
2:  Residue( ALA:53  nAtoms=5 )
3:  Residue( ALA:65  nAtoms=5 )
4:  Residue( ALA:85  nAtoms=5 )
)

gives the first 5 residues named "ALA".

Searching for residues
----------------------

You can also search for residues, using their name (``resname``),
their number (``resnum``) and/or their index in their parent
molecule (``residx``).

>>> print(mol.residues("resnum 5"))
Selector<SireMol::Residue>( size=2
0:  Residue( GLU:5   nAtoms=9 )
1:  Residue( GLU:5   nAtoms=9 )
)

.. note::

    There are two residues with number 5 as there are multiple chains
    in this protein. Note also how the residue's name (GLU) and
    number (5) are printed in its output.

You can use the residue search string in a molecular container's index
operator too!

>>> print(mol["resnum 5"])
Selector<SireMol::Residue>( size=2
0:  Residue( GLU:5   nAtoms=9 )
1:  Residue( GLU:5   nAtoms=9 )
)

and you can combine it with atom identifiers, e.g.

>>> print(mol["resname ALA and atomname CA"])
Selector<SireMol::Atom>( size=155
0:  Atom( CA:65   [ -54.77,   13.35,   37.26] )
1:  Atom( CA:117  [ -62.33,   13.58,   32.15] )
2:  Atom( CA:204  [ -45.04,    6.02,   36.66] )
3:  Atom( CA:306  [ -47.63,   28.39,   36.61] )
4:  Atom( CA:352  [ -34.57,   20.94,   29.60] )
...
150:  Atom( CA:10774 [  -4.40,    7.58,   14.84] )
151:  Atom( CA:10816 [  -1.17,    9.47,   25.09] )
152:  Atom( CA:10886 [   9.70,  -11.41,   19.28] )
153:  Atom( CA:11247 [  14.11,    2.16,   14.69] )
154:  Atom( CA:11624 [  22.43,   -6.30,   32.21] )
)

You can also search for multiple residue names or numbers.

>>> print(mol["resname ALA, ARG"])
Selector<SireMol::Residue>( size=255
0:  Residue( ALA:23  nAtoms=5 )
1:  Residue( ALA:30  nAtoms=5 )
2:  Residue( ALA:53  nAtoms=5 )
3:  Residue( ARG:61  nAtoms=11 )
4:  Residue( ALA:65  nAtoms=5 )
...
250:  Residue( ARG:652 nAtoms=11 )
251:  Residue( ARG:657 nAtoms=11 )
252:  Residue( ARG:680 nAtoms=11 )
253:  Residue( ARG:685 nAtoms=11 )
254:  Residue( ALA:691 nAtoms=5 )
)

>>> print(mol["resnum 5, 7, 9"])
Selector<SireMol::Residue>( size=10
0:  Residue( VAL:7   nAtoms=7 )
1:  Residue( LYS:9   nAtoms=9 )
2:  Residue( GLU:5   nAtoms=9 )
3:  Residue( VAL:7   nAtoms=7 )
4:  Residue( GLU:9   nAtoms=9 )
5:  Residue( VAL:7   nAtoms=7 )
6:  Residue( LYS:9   nAtoms=9 )
7:  Residue( GLU:5   nAtoms=9 )
8:  Residue( VAL:7   nAtoms=7 )
9:  Residue( GLU:9   nAtoms=9 )
)

Finding the atoms in a residue
------------------------------

Because both :class:`~Sire.Mol.Residue` and :class:`~Sire.Mol.Selector_Residue_`
are molecular containers, they also have their own
:func:`~Sire.Mol.Residue.atom` and :func:`~Sire.Mol.Residue.atoms` functions,
which behave as you would expect.

>>> print(mol["resname ALA"].atoms("CA"))
Selector<SireMol::Atom>( size=155
0:  Atom( CA:65   [ -54.77,   13.35,   37.26] )
1:  Atom( CA:117  [ -62.33,   13.58,   32.15] )
2:  Atom( CA:204  [ -45.04,    6.02,   36.66] )
3:  Atom( CA:306  [ -47.63,   28.39,   36.61] )
4:  Atom( CA:352  [ -34.57,   20.94,   29.60] )
...
150:  Atom( CA:10774 [  -4.40,    7.58,   14.84] )
151:  Atom( CA:10816 [  -1.17,    9.47,   25.09] )
152:  Atom( CA:10886 [   9.70,  -11.41,   19.28] )
153:  Atom( CA:11247 [  14.11,    2.16,   14.69] )
154:  Atom( CA:11624 [  22.43,   -6.30,   32.21] )
)

You can get all of the atoms in a residue by calling the
:func:`~Sire.Mol.Residue.atoms` function without any arguments.

>>> Selector<SireMol::Atom>( size=8
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
2:  Atom( C:3     [ -56.06,    9.95,   42.55] )
3:  Atom( O:4     [ -57.04,    9.73,   41.82] )
4:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
5:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
6:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
7:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
)

Uniquely identifying a residue
------------------------------

You uniquely identify a residue in a molecule using its residue index
(``residx``). You can get the index of a residue in a molecule by
calling its :func:`~Sire.Mol.Residue.index` function.

>>> print(mol.residue(0).index())
ResIdx(0)

Iterating over residues
-----------------------

The :class:`~Sire.Mol.Selector_Residue_` class is iterable, meaning that
it can be used in loops.

>>> for res in mol["resname ALA and resnum < 30"]:
...     print(res)
Residue( ALA:23  nAtoms=5 )
Residue( ALA:16  nAtoms=5 )
Residue( ALA:21  nAtoms=5 )
Residue( ALA:23  nAtoms=5 )
Residue( ALA:16  nAtoms=5 )

This is particulary useful when combined with looping over the
atoms in the residues.

>>> for res in mol["residx < 3"]:
...     for atom in res["atomname C, CA"]:
...         print(res, atom)
Residue( ILE:6   nAtoms=8 ) Atom( CA:2    [ -55.43,   11.35,   42.54] )
Residue( ILE:6   nAtoms=8 ) Atom( C:3     [ -56.06,    9.95,   42.55] )
Residue( VAL:7   nAtoms=7 ) Atom( CA:10   [ -56.02,    7.64,   43.47] )
Residue( VAL:7   nAtoms=7 ) Atom( C:11    [ -56.14,    7.05,   42.06] )
Residue( LEU:8   nAtoms=8 ) Atom( CA:17   [ -54.99,    6.39,   39.98] )
Residue( LEU:8   nAtoms=8 ) Atom( C:18    [ -54.61,    4.90,   40.03] )
