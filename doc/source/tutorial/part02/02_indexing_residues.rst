=================
Indexing Residues
=================

Residues are collections of atoms. They typically represent an amino
acid residue in a protein. Residues are implemented via the
:class:`~sire.mol.Residue` class, which itself is a molecular container
for :class:`~sire.mol.Atom` objects. An atom can only belong to one
residue at a time (and they don't need to be assigned to a residue).

You can access residues in a molecule container using the
:func:`~sire.mol.Residue.residue` and :func:`~sire.mol.Residue.residues`
functions, which are available on all of the molecular container types.

>>> print(mol.residue(0))
Residue( ILE:6   num_atoms=8 )

gives the molecule at index 0, while

>>> print(mol.residues("ALA"))
Selector<SireMol::Residue>( size=155
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
...
150:  Residue( ALA:578 num_atoms=5 )
151:  Residue( ALA:584 num_atoms=5 )
152:  Residue( ALA:593 num_atoms=5 )
153:  Residue( ALA:646 num_atoms=5 )
154:  Residue( ALA:691 num_atoms=5 )
)

returns all residues that are named "ALA".

The :func:`~sire.mol.Residue.residue` will raise a KeyError if more than
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
0:  Residue( ILE:6   num_atoms=8 )
1:  Residue( VAL:7   num_atoms=7 )
2:  Residue( LEU:8   num_atoms=8 )
3:  Residue( LYS:9   num_atoms=9 )
4:  Residue( SER:10  num_atoms=6 )
5:  Residue( SER:11  num_atoms=6 )
6:  Residue( ASP:12  num_atoms=8 )
7:  Residue( GLY:13  num_atoms=4 )
8:  Residue( VAL:22  num_atoms=7 )
9:  Residue( ALA:23  num_atoms=5 )
)

The result, a :class:`~sire.mol.Selector_Residue_` is also a molecular
container, and can be used like :class:`~sire.mol.Selector_Atom_`.

>>> print(mol.residues("ALA")[0:5])
Selector<SireMol::Residue>( size=5
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
)

gives the first 5 residues named "ALA".

Searching for residues
----------------------

You can also search for residues, using their name (``resname``),
their number (``resnum``) and/or their index in their parent
molecule (``residx``).

>>> print(mol.residues("resnum 5"))
Selector<SireMol::Residue>( size=2
0:  Residue( GLU:5   num_atoms=9 )
1:  Residue( GLU:5   num_atoms=9 )
)

.. note::

    There are two residues with number 5 as there are multiple chains
    in this protein. Note also how the residue's name (GLU) and
    number (5) are printed in its output.

You can use the residue search string in a molecular container's index
operator too!

>>> print(mol["resnum 5"])
Selector<SireMol::Residue>( size=2
0:  Residue( GLU:5   num_atoms=9 )
1:  Residue( GLU:5   num_atoms=9 )
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
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ARG:61  num_atoms=11 )
4:  Residue( ALA:65  num_atoms=5 )
...
250:  Residue( ARG:652 num_atoms=11 )
251:  Residue( ARG:657 num_atoms=11 )
252:  Residue( ARG:680 num_atoms=11 )
253:  Residue( ARG:685 num_atoms=11 )
254:  Residue( ALA:691 num_atoms=5 )
)

>>> print(mol["resnum 5, 7, 9"])
Selector<SireMol::Residue>( size=10
0:  Residue( VAL:7   num_atoms=7 )
1:  Residue( LYS:9   num_atoms=9 )
2:  Residue( GLU:5   num_atoms=9 )
3:  Residue( VAL:7   num_atoms=7 )
4:  Residue( GLU:9   num_atoms=9 )
5:  Residue( VAL:7   num_atoms=7 )
6:  Residue( LYS:9   num_atoms=9 )
7:  Residue( GLU:5   num_atoms=9 )
8:  Residue( VAL:7   num_atoms=7 )
9:  Residue( GLU:9   num_atoms=9 )
)

>>> print(mol["resnum 201:205"])
Selector<SireMol::Residue>( size=9
0:  Residue( LEU:201 num_atoms=8 )
1:  Residue( ARG:202 num_atoms=11 )
2:  Residue( GLU:203 num_atoms=9 )
3:  Residue( LEU:204 num_atoms=8 )
4:  Residue( LEU:201 num_atoms=8 )
5:  Residue( ARG:202 num_atoms=11 )
6:  Residue( GLU:203 num_atoms=9 )
7:  Residue( LEU:204 num_atoms=8 )
8:  Residue( PEG:201 num_atoms=7 )
)

Wildcard (glob) searching is also supported for residue names.

>>> print(mol["resname /ala/i"])
Selector<SireMol::Residue>( size=155
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
...
150:  Residue( ALA:578 num_atoms=5 )
151:  Residue( ALA:584 num_atoms=5 )
152:  Residue( ALA:593 num_atoms=5 )
153:  Residue( ALA:646 num_atoms=5 )
154:  Residue( ALA:691 num_atoms=5 )
)

>>> print(mol["resname /HI?/"])
Selector<SireMol::Residue>( size=42
0:  Residue( HIS:62  num_atoms=10 )
1:  Residue( HIS:27  num_atoms=10 )
2:  Residue( HIS:39  num_atoms=10 )
3:  Residue( HIS:75  num_atoms=10 )
4:  Residue( HIS:84  num_atoms=10 )
...
37:  Residue( HIS:638 num_atoms=10 )
38:  Residue( HIS:639 num_atoms=10 )
39:  Residue( HIS:662 num_atoms=10 )
40:  Residue( HIS:666 num_atoms=10 )
41:  Residue( HIS:668 num_atoms=10 )
)

This last search is particularly useful for proteins, as it is common
for histidine residues to have different names depending on protonation
state (e.g. "HIS", "HIP", "HIE" or "HID").

Finding the atoms in a residue
------------------------------

Because both :class:`~sire.mol.Residue` and :class:`~sire.mol.Selector_Residue_`
are molecular containers, they also have their own
:func:`~sire.mol.Residue.atom` and :func:`~sire.mol.Residue.atoms` functions,
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
:func:`~sire.mol.Residue.atoms` function without any arguments.

>>> mol["residx 0"].atoms()
Selector<SireMol::Atom>( size=8
0:  Atom( N:1     [ -54.07,   11.27,   41.93] )
1:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
2:  Atom( C:3     [ -56.06,    9.95,   42.55] )
3:  Atom( O:4     [ -57.04,    9.73,   41.82] )
4:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
5:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
6:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
7:  Atom( CD1:8   [ -55.42,   14.31,   43.09] )
)

Another route is to use the ``atoms in`` phrase in a search string, e.g.

>>> print(mol["atoms in resname ALA"])
Selector<SireMol::Atom>( size=775
0:  Atom( N:64    [ -54.11,   14.36,   38.13] )
1:  Atom( CA:65   [ -54.77,   13.35,   37.26] )
2:  Atom( C:66    [ -55.92,   14.01,   36.49] )
3:  Atom( O:67    [ -57.09,   13.65,   36.74] )
4:  Atom( CB:68   [ -55.25,   12.19,   38.09] )
...
770:  Atom( N:11623 [  22.09,   -7.64,   32.65] )
771:  Atom( CA:11624 [  22.43,   -6.30,   32.21] )
772:  Atom( C:11625 [  23.84,   -6.28,   31.63] )
773:  Atom( O:11626 [  24.72,   -7.01,   32.08] )
774:  Atom( CB:11627 [  22.32,   -5.30,   33.36] )
)

This has returned all of the atoms in residues that are called "ALA".

You can get the residues that match atoms using ``residues with``, e.g.

>>> print(mol["residues with atomname CA"])
Selector<SireMol::Residue>( size=1494
0:  Residue( ILE:6   num_atoms=8 )
1:  Residue( VAL:7   num_atoms=7 )
2:  Residue( LEU:8   num_atoms=8 )
3:  Residue( LYS:9   num_atoms=9 )
4:  Residue( SER:10  num_atoms=6 )
...
1489:  Residue( ALA:691 num_atoms=5 )
1490:  Residue( PRO:692 num_atoms=7 )
1491:  Residue( GLU:693 num_atoms=9 )
1492:  Residue( ASN:694 num_atoms=8 )
1493:  Residue( ASP:695 num_atoms=8 )
)

This has returned all of the residues that contain an atom called "CA".

Another way to do this would be to call the :func:`~sire.mol.Selected_Atom_.residues`
function on the molecular container, e.g.

>>> print(mol["CA"].residues())
Selector<SireMol::Residue>( size=1494
0:  Residue( ILE:6   num_atoms=8 )
1:  Residue( VAL:7   num_atoms=7 )
2:  Residue( LEU:8   num_atoms=8 )
3:  Residue( LYS:9   num_atoms=9 )
4:  Residue( SER:10  num_atoms=6 )
...
1489:  Residue( ALA:691 num_atoms=5 )
1490:  Residue( PRO:692 num_atoms=7 )
1491:  Residue( GLU:693 num_atoms=9 )
1492:  Residue( ASN:694 num_atoms=8 )
1493:  Residue( ASP:695 num_atoms=8 )
)

Uniquely identifying a residue
------------------------------

You uniquely identify a residue in a molecule using its residue index
(``residx``). You can get the index of a residue in a molecule by
calling its :func:`~sire.mol.Residue.index` function.

>>> print(mol.residue(0).index())
ResIdx(0)

.. warning::

    Be careful indexing by residue index. This is the index of the residue
    that uniquely identifies it within its parent molecule. It is not the
    index of the residue in an arbitrary molecular container.

Residue identifying types
-------------------------

Another way to index residues is to use the residue indexing types, i.e.
:class:`~sire.mol.ResIdx`, :class:`~sire.mol.ResName` and
:class:`~sire.mol.ResNum`. The easiest way to create these is
by using the function :func:`sire.resid`.

>>> print(mol[sr.resid("ALA")])
Selector<SireMol::Residue>( size=155
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:30  num_atoms=5 )
2:  Residue( ALA:53  num_atoms=5 )
3:  Residue( ALA:65  num_atoms=5 )
4:  Residue( ALA:85  num_atoms=5 )
...
150:  Residue( ALA:578 num_atoms=5 )
151:  Residue( ALA:584 num_atoms=5 )
152:  Residue( ALA:593 num_atoms=5 )
153:  Residue( ALA:646 num_atoms=5 )
154:  Residue( ALA:691 num_atoms=5 )
)

This returns the residues called "ALA", as ``sr.resid("ALA")`` has created
an :class:`~sire.mol.ResName` object.

>>> print(sr.resid("ALA"))
ResName('ALA')

This function will create an :class:`~sire.mol.ResNum` if it is passed
an integer, e.g.

>>> print(sr.resid(5))
ResNum(5)
>>> print(mol[sr.resid(5)])
Selector<SireMol::Residue>( size=2
0:  Residue( GLU:5   num_atoms=9 )
1:  Residue( GLU:5   num_atoms=9 )
)

You can set both a name and a number by passing in two arguments, e.g.

>>> print(mol[sr.resid("ALA", 23)])
Selector<SireMol::Residue>( size=2
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:23  num_atoms=5 )
)
>>> print(mol[sr.resid(name="ALA", num=23)])
Selector<SireMol::Residue>( size=2
0:  Residue( ALA:23  num_atoms=5 )
1:  Residue( ALA:23  num_atoms=5 )
)

Iterating over residues
-----------------------

The :class:`~sire.mol.Selector_Residue_` class is iterable, meaning that
it can be used in loops.

>>> for res in mol["resname ALA and resnum < 30"]:
...     print(res)
Residue( ALA:23  num_atoms=5 )
Residue( ALA:16  num_atoms=5 )
Residue( ALA:21  num_atoms=5 )
Residue( ALA:23  num_atoms=5 )
Residue( ALA:16  num_atoms=5 )

This is particulary useful when combined with looping over the
atoms in the residues.

>>> for res in mol["residx < 3"]:
...     for atom in res["atomname C, CA"]:
...         print(res, atom)
Residue( ILE:6   num_atoms=8 ) Atom( CA:2    [ -55.43,   11.35,   42.54] )
Residue( ILE:6   num_atoms=8 ) Atom( C:3     [ -56.06,    9.95,   42.55] )
Residue( VAL:7   num_atoms=7 ) Atom( CA:10   [ -56.02,    7.64,   43.47] )
Residue( VAL:7   num_atoms=7 ) Atom( C:11    [ -56.14,    7.05,   42.06] )
Residue( LEU:8   num_atoms=8 ) Atom( CA:17   [ -54.99,    6.39,   39.98] )
Residue( LEU:8   num_atoms=8 ) Atom( C:18    [ -54.61,    4.90,   40.03] )

Counting residues
-----------------

Similar to how you did for atom, you can find the set of residue names
via

>>> print(set(mol.residues().names()))
{ResName('ALA'),
 ResName('ARG'),
 ResName('ASN'),
 ResName('ASP'),
 ResName('CIT'),
 ResName('CYS'),
 ResName('GLN'),
 ResName('GLU'),
 ResName('GLY'),
 ResName('HIS'),
 ResName('HOH'),
 ResName('ILE'),
 ResName('LEU'),
 ResName('LYS'),
 ResName('MET'),
 ResName('PEG'),
 ResName('PHE'),
 ResName('PRO'),
 ResName('SER'),
 ResName('THR'),
 ResName('TRP'),
 ResName('TYR'),
 ResName('VAL')}

And you can count how many of each residue using;

>>> for name in set(mol.residues().names()):
...     print(name, len(mol.residues(name)))
ResName('VAL') 74
ResName('ILE') 64
ResName('GLN') 32
ResName('PRO') 90
ResName('GLU') 107
ResName('TRP') 24
ResName('GLY') 68
ResName('CYS') 48
ResName('HOH') 18
ResName('CIT') 2
ResName('ARG') 100
ResName('MET') 20
ResName('SER') 102
ResName('PHE') 64
ResName('ASN') 38
ResName('THR') 88
ResName('ASP') 84
ResName('LYS') 46
ResName('TYR') 22
ResName('HIS') 42
ResName('PEG') 4
ResName('ALA') 155
ResName('LEU') 226

This can be a convenient way of finding the residue names of different
ligands or cofactors that are bound to the molecule.

You could do a similar thing for residue numbers, e.g.

>>> for number in set(mol.residues().numbers()):
...     print(number, len(mol.residues(number)))
ResNum(5) 2
ResNum(6) 4
ResNum(7) 4
ResNum(8) 4
ResNum(9) 4
ResNum(10) 4
ResNum(11) 4
ResNum(12) 4
ResNum(13) 4
ResNum(14) 2
...
