==================
Indexing Molecules
==================

Molecules are collections of atoms that (at least notionally) should all
be bonded to each other and represent the concept of a molecule in a
chemical system.

A molecule is a molecular container for atoms, residues, chains and
segments, and is implemented via the :class:`sire.mol.Molecule` class.
At a minimum, a :class:`sire.mol.Molecule` contains at least one
:class:`sire.mol.Atom`. There
is no requirement for this atom to be placed into a residue, chain
or segment, and so it is possible for a molecule to have zero
residues, chains or segments.

There are several classes in ``sire`` that can act as molecular containers
for molecule objects. You have already encountered one, which is the
:class:`~sire.system.System` class, into which the molecules were first
loaded.

For example, load up the ``aladip`` system;

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> print(mols)
System( name=ACE num_molecules=631 num_residues=633 num_atoms=1912 )

This has loaded 631 molecules. You can access molecules by their index
in the system by simple index;

>>> print(mols[0])
Molecule( ACE:2   num_atoms=22 num_residues=3 )

or by their name

>>> print(mols["ACE"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )

You can also index via a range

>>> print(mols[0:10])
SelectorMol( size=10
0: Molecule( ACE:2   num_atoms=22 num_residues=3 )
1: Molecule( WAT:3   num_atoms=3 num_residues=1 )
2: Molecule( WAT:4   num_atoms=3 num_residues=1 )
3: Molecule( WAT:5   num_atoms=3 num_residues=1 )
4: Molecule( WAT:6   num_atoms=3 num_residues=1 )
5: Molecule( WAT:7   num_atoms=3 num_residues=1 )
6: Molecule( WAT:8   num_atoms=3 num_residues=1 )
7: Molecule( WAT:9   num_atoms=3 num_residues=1 )
8: Molecule( WAT:10  num_atoms=3 num_residues=1 )
9: Molecule( WAT:11  num_atoms=3 num_residues=1 )
)

or via a list of indicies

>>> print(mols[ [0, 2, 4, 6] ])
SelectorMol( size=4
0: Molecule( ACE:2   num_atoms=22 num_residues=3 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:6   num_atoms=3 num_residues=1 )
3: Molecule( WAT:8   num_atoms=3 num_residues=1 )
)

If multiple molecules are returned, then the result is held in a
:class:`~sire.mol.SelectorMol` class. This is also a molecular container,
and can be used in the same way as all other molecular containers. For example,

>>> print(mols["WAT"][0])
Molecule( WAT:3   num_atoms=3 num_residues=1 )

gives the first molecule called ``WAT``.

You can achieve the same result using the :func:`~sire.system.System.molecule`
and :func:`~sire.system.System.molecules` functions.

>>> print(mols.molecule(0))
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols.molecules(range(0, 10)))
SelectorMol( size=10
0: Molecule( ACE:2   num_atoms=22 num_residues=3 )
1: Molecule( WAT:3   num_atoms=3 num_residues=1 )
2: Molecule( WAT:4   num_atoms=3 num_residues=1 )
3: Molecule( WAT:5   num_atoms=3 num_residues=1 )
4: Molecule( WAT:6   num_atoms=3 num_residues=1 )
5: Molecule( WAT:7   num_atoms=3 num_residues=1 )
6: Molecule( WAT:8   num_atoms=3 num_residues=1 )
7: Molecule( WAT:9   num_atoms=3 num_residues=1 )
8: Molecule( WAT:10  num_atoms=3 num_residues=1 )
9: Molecule( WAT:11  num_atoms=3 num_residues=1 )
)
>>> print(mols.molecules("WAT").molecule(0))
Molecule( WAT:3   num_atoms=3 num_residues=1 )

Search for molecules
--------------------

You can search for a molecule in a molecule container via their
name (``molname``) or number (``molnum``).

>>> print(mols["molname ACE"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols["molnum 2"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols["molidx 0"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols["molname WAT"])
SelectorMol( size=630
0: Molecule( WAT:3   num_atoms=3 num_residues=1 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:5   num_atoms=3 num_residues=1 )
3: Molecule( WAT:6   num_atoms=3 num_residues=1 )
4: Molecule( WAT:7   num_atoms=3 num_residues=1 )
...
625: Molecule( WAT:628 num_atoms=3 num_residues=1 )
626: Molecule( WAT:629 num_atoms=3 num_residues=1 )
627: Molecule( WAT:630 num_atoms=3 num_residues=1 )
628: Molecule( WAT:631 num_atoms=3 num_residues=1 )
629: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

.. note::

   Note how the name and number of the molecule are printed out,
   e.g. ``ACE:2`` is for the molecule called ``ACE`` with number ``2``.

You can also search by the index in the molecular container (``molidx``),
e.g.

>>> print(mols["molidx > 50"])
SelectorMol( size=580
0: Molecule( WAT:53  num_atoms=3 num_residues=1 )
1: Molecule( WAT:54  num_atoms=3 num_residues=1 )
2: Molecule( WAT:55  num_atoms=3 num_residues=1 )
3: Molecule( WAT:56  num_atoms=3 num_residues=1 )
4: Molecule( WAT:57  num_atoms=3 num_residues=1 )
...
575: Molecule( WAT:628 num_atoms=3 num_residues=1 )
576: Molecule( WAT:629 num_atoms=3 num_residues=1 )
577: Molecule( WAT:630 num_atoms=3 num_residues=1 )
578: Molecule( WAT:631 num_atoms=3 num_residues=1 )
579: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

.. note::

   The ``molidx`` is the index of the molecule in its parent container.
   This will depend on which container the search is performed with.
   This is different to ``atomidx``, ``residx`` etc, which are unique
   identifying indicies of the atom, residue etc in their parent molecule.

You can combine the search with other identifiers, e.g.

>>> print(mols["molname WAT and molnum 3"])
Molecule( WAT:3   num_atoms=3 num_residues=1 )

and can search for multiple names or numbers

>>> print(mols["molname ACE, WAT"])
SelectorMol( size=631
0: Molecule( ACE:2   num_atoms=22 num_residues=3 )
1: Molecule( WAT:3   num_atoms=3 num_residues=1 )
2: Molecule( WAT:4   num_atoms=3 num_residues=1 )
3: Molecule( WAT:5   num_atoms=3 num_residues=1 )
4: Molecule( WAT:6   num_atoms=3 num_residues=1 )
...
626: Molecule( WAT:628 num_atoms=3 num_residues=1 )
627: Molecule( WAT:629 num_atoms=3 num_residues=1 )
628: Molecule( WAT:630 num_atoms=3 num_residues=1 )
629: Molecule( WAT:631 num_atoms=3 num_residues=1 )
630: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)
>>> print(mols["molnum 5:11, 30, 40"])
SelectorMol( size=8
0: Molecule( WAT:5   num_atoms=3 num_residues=1 )
1: Molecule( WAT:6   num_atoms=3 num_residues=1 )
2: Molecule( WAT:7   num_atoms=3 num_residues=1 )
3: Molecule( WAT:8   num_atoms=3 num_residues=1 )
4: Molecule( WAT:9   num_atoms=3 num_residues=1 )
5: Molecule( WAT:10  num_atoms=3 num_residues=1 )
6: Molecule( WAT:30  num_atoms=3 num_residues=1 )
7: Molecule( WAT:40  num_atoms=3 num_residues=1 )
)

Wildcard (glob) searching is also supported for molecule names.

>>> print(mols["molname /A*/"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols["molname /wat/i"])
SelectorMol( size=630
0: Molecule( WAT:3   num_atoms=3 num_residues=1 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:5   num_atoms=3 num_residues=1 )
3: Molecule( WAT:6   num_atoms=3 num_residues=1 )
4: Molecule( WAT:7   num_atoms=3 num_residues=1 )
...
625: Molecule( WAT:628 num_atoms=3 num_residues=1 )
626: Molecule( WAT:629 num_atoms=3 num_residues=1 )
627: Molecule( WAT:630 num_atoms=3 num_residues=1 )
628: Molecule( WAT:631 num_atoms=3 num_residues=1 )
629: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

Finding atoms, residues, chains and segments in a molecule
----------------------------------------------------------

Because both :class:`~sire.system.System` and :class:`~sire.mol.SelectorMol`
are molecular containers, they both have their own
:func:`~sire.mol.SelectorMol.atom`, :func:`~sire.mol.SelectorMol.atoms`,
:func:`~sire.mol.SelectorMol.residue`, :func:`~sire.mol.SelectorMol.residues`,
etc functions for accessing atoms, residues, chains or segments across
multiple molecules.

For example, you can get all of the atoms that are oxygens using

>>> print(mols.atoms("element O"))
SireMol::SelectorM<SireMol::Atom>( size=632
0: MolNum(2) Atom( O:6     [  19.19,    5.44,   14.76] )
1: MolNum(2) Atom( O:16    [  14.94,    3.17,   15.88] )
2: MolNum(3) Atom( O:23    [  25.64,    8.50,   22.42] )
3: MolNum(4) Atom( O:26    [  22.83,    8.93,    4.14] )
4: MolNum(5) Atom( O:29    [   9.96,   24.84,   23.52] )
...
627: MolNum(628) Atom( O:1898  [  22.40,   11.84,   10.08] )
628: MolNum(629) Atom( O:1901  [   0.63,    8.69,   19.94] )
629: MolNum(630) Atom( O:1904  [  18.69,   22.12,    9.35] )
630: MolNum(631) Atom( O:1907  [  21.65,    7.88,    9.79] )
631: MolNum(632) Atom( O:1910  [   9.25,   13.73,    2.29] )
)

or, even easier,

>>> print(mols["element O"])
SireMol::SelectorM<SireMol::Atom>( size=632
0: MolNum(2) Atom( O:6     [  19.19,    5.44,   14.76] )
1: MolNum(2) Atom( O:16    [  14.94,    3.17,   15.88] )
2: MolNum(3) Atom( O:23    [  25.64,    8.50,   22.42] )
3: MolNum(4) Atom( O:26    [  22.83,    8.93,    4.14] )
4: MolNum(5) Atom( O:29    [   9.96,   24.84,   23.52] )
...
627: MolNum(628) Atom( O:1898  [  22.40,   11.84,   10.08] )
628: MolNum(629) Atom( O:1901  [   0.63,    8.69,   19.94] )
629: MolNum(630) Atom( O:1904  [  18.69,   22.12,    9.35] )
630: MolNum(631) Atom( O:1907  [  21.65,    7.88,    9.79] )
631: MolNum(632) Atom( O:1910  [   9.25,   13.73,    2.29] )
)

The result is a :class:`~sire.mol.SelectorM_Atom_`. This is a multi-molecule
version of :class:`~sire.mol.Selector_Atom_`, which behaves in an identical
way.

This works for residues too!

>>> print(mols.residues("WAT"))
SireMol::SelectorM<SireMol::Residue>( size=630
0: MolNum(3) Residue( WAT:4   num_atoms=3 )
1: MolNum(4) Residue( WAT:5   num_atoms=3 )
2: MolNum(5) Residue( WAT:6   num_atoms=3 )
3: MolNum(6) Residue( WAT:7   num_atoms=3 )
4: MolNum(7) Residue( WAT:8   num_atoms=3 )
...
625: MolNum(628) Residue( WAT:629 num_atoms=3 )
626: MolNum(629) Residue( WAT:630 num_atoms=3 )
627: MolNum(630) Residue( WAT:631 num_atoms=3 )
628: MolNum(631) Residue( WAT:632 num_atoms=3 )
629: MolNum(632) Residue( WAT:633 num_atoms=3 )
)
>>> print(mols["resname WAT"])
SelectorMol( size=630
0: Molecule( WAT:3   num_atoms=3 num_residues=1 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:5   num_atoms=3 num_residues=1 )
3: Molecule( WAT:6   num_atoms=3 num_residues=1 )
4: Molecule( WAT:7   num_atoms=3 num_residues=1 )
...
625: Molecule( WAT:628 num_atoms=3 num_residues=1 )
626: Molecule( WAT:629 num_atoms=3 num_residues=1 )
627: Molecule( WAT:630 num_atoms=3 num_residues=1 )
628: Molecule( WAT:631 num_atoms=3 num_residues=1 )
629: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

.. note::

   Note how the index operator will convert the result of the search to
   the largest possible container (in this case, ``Molecule``), while
   the ``.residues`` function will always return the result as a
   set of ``Residue`` objects.

Another route is to use the ``atoms in``, ``residues in``, ``chains in``, or
``segments in`` phrases in the search, e.g.

>>> print(mols["atoms in molname ACE"])
Selector<SireMol::Atom>( size=22
0:  Atom( HH31:1  [  18.45,    3.49,   12.44] )
1:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
2:  Atom( HH32:3  [  20.05,    3.63,   13.29] )
3:  Atom( HH33:4  [  18.80,    2.43,   13.73] )
4:  Atom( C:5     [  18.48,    4.55,   14.35] )
...
17:  Atom( H:18    [  15.34,    5.45,   17.96] )
18:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
19:  Atom( HH31:20 [  14.35,    3.41,   19.15] )
20:  Atom( HH32:21 [  13.19,    4.59,   18.94] )
21:  Atom( HH33:22 [  13.21,    3.33,   17.69] )
)

.. note::

   Note how the index operator will return a single molecule ``Selector_Atom_``
   when only a single molecule matches the search.

>>> print(mols["residues in molnum 10:21"])
SelectorMol( size=11
0: Molecule( WAT:10  num_atoms=3 num_residues=1 )
1: Molecule( WAT:11  num_atoms=3 num_residues=1 )
2: Molecule( WAT:12  num_atoms=3 num_residues=1 )
3: Molecule( WAT:13  num_atoms=3 num_residues=1 )
4: Molecule( WAT:14  num_atoms=3 num_residues=1 )
...
6: Molecule( WAT:16  num_atoms=3 num_residues=1 )
7: Molecule( WAT:17  num_atoms=3 num_residues=1 )
8: Molecule( WAT:18  num_atoms=3 num_residues=1 )
9: Molecule( WAT:19  num_atoms=3 num_residues=1 )
10: Molecule( WAT:20  num_atoms=3 num_residues=1 )
)

You can also go back to the containing molecule using ``molecules with``,

>>> print(mols["molecules with element C"])
Molecule( ACE:2   num_atoms=22 num_residues=3 )
>>> print(mols["molecules with resname /wat/i"])
SelectorMol( size=630
0: Molecule( WAT:3   num_atoms=3 num_residues=1 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:5   num_atoms=3 num_residues=1 )
3: Molecule( WAT:6   num_atoms=3 num_residues=1 )
4: Molecule( WAT:7   num_atoms=3 num_residues=1 )
...
625: Molecule( WAT:628 num_atoms=3 num_residues=1 )
626: Molecule( WAT:629 num_atoms=3 num_residues=1 )
627: Molecule( WAT:630 num_atoms=3 num_residues=1 )
628: Molecule( WAT:631 num_atoms=3 num_residues=1 )
629: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

This last search is particularly useful if you are looking for a protein
or for a ligand molecule in a system. In these cases, you could look for
``molecules with resname /ala/i`` or ``molecules with resname /lig/i``.

Uniquely identifying a molecule
-------------------------------

Molecules are uniquely identified by their molecule number. This is a number
that ``sire`` assigns to molecules when they are loaded. ``sire`` ensures
that each molecule that it loads will have its own unique number. You have
no control over the number, and should not assume that the number will be
the same every time you run your script. The number is the count of the
molecules that have been loaded during a ``sire`` session, e.g. molecule
number 2 refers to the second molecule that has been loaded. You can assume
that the first molecule loaded from a file will have the smallest molecule
number.

You can get all of the numbers of the molecules in a container using

>>> print(mols.numbers())
[MolNum(2), MolNum(3), MolNum(4), MolNum(5)... MolNum(632)]

.. note::

   You can use ``mols.names()`` to get all of the molecule names.

``sire`` uses the molecule number as it is the only identifier that can
be guaranteed to be unique in the program. Molecules can be moved and copied
into multiple containers, and its index will depend on its container.
Molecule names are assigned from the file, and multiple molecules can
have the same name.

Molecule identifying types
--------------------------

Another way to index molecules is to use the molecule identifying types, i.e.
:class:`~sire.mol.MolName`, :class:`~sire.mol.MolNum` and
:class:`~sire.mol.MolIdx`. The easiest way to create these is via the
function :func:`sire.molid`.

>>> print(mols[sr.molid("ACE")])
Molecule( ACE:2   num_atoms=22 num_residues=3 )

This returns the molecule called "ACE", as ``sr.molid("ACE")`` has created
an :class:`~sire.mol.MolName` object.

>>> print(sr.molid("ACE"))
MolName('ACE')

This function will create an :class:`~sire.mol.MolNum` if it is passed
an integer, e.g.

>>> print(sr.molid(5))
MolNum(5)
>>> print(mols[sr.molid(5)])
Molecule( WAT:5   num_atoms=3 num_residues=1 )

You can set both a name and a number by passing in two arguments, e.g.

>>> print(mols[sr.molid("WAT", 5)])
Molecule( WAT:5   num_atoms=3 num_residues=1 )

Iterating over molecules
------------------------

The :class:`~sire.mol.SelectorMol` and :class:`~sire.system.System` classes
are iterable, meaning that they can be used in loops.

>>> for mol in mols["molname WAT and molnum < 10"]:
...     print(mol)
Molecule( WAT:3   num_atoms=3 num_residues=1 )
Molecule( WAT:4   num_atoms=3 num_residues=1 )
Molecule( WAT:5   num_atoms=3 num_residues=1 )
Molecule( WAT:6   num_atoms=3 num_residues=1 )
Molecule( WAT:7   num_atoms=3 num_residues=1 )
Molecule( WAT:8   num_atoms=3 num_residues=1 )
Molecule( WAT:9   num_atoms=3 num_residues=1 )

This is particulary useful when combined with looping over the
atoms or residues in the molecules.

>>> for mol in mols["molidx < 3"]:
...     for atom in mol["element O"]:
...         print(mol, atom)
Molecule( ACE:2   num_atoms=22 num_residues=3 ) Atom( O:6     [  19.19,    5.44,   14.76] )
Molecule( ACE:2   num_atoms=22 num_residues=3 ) Atom( O:16    [  14.94,    3.17,   15.88] )
Molecule( WAT:3   num_atoms=3 num_residues=1 ) Atom( O:23    [  25.64,    8.50,   22.42] )
Molecule( WAT:4   num_atoms=3 num_residues=1 ) Atom( O:26    [  22.83,    8.93,    4.14] )

Smart search terms
------------------

There are a few smart search terms that can help find molecules.

The most useful is perhaps ``water``. This searches for molecules that
contain one oxygen atom, two hydrogen atoms and any number of null
element (dummy) atoms.

>>> print(mols["water"])
SelectorMol( size=630
0: Molecule( WAT:3   num_atoms=3 num_residues=1 )
1: Molecule( WAT:4   num_atoms=3 num_residues=1 )
2: Molecule( WAT:5   num_atoms=3 num_residues=1 )
3: Molecule( WAT:6   num_atoms=3 num_residues=1 )
4: Molecule( WAT:7   num_atoms=3 num_residues=1 )
...
625: Molecule( WAT:628 num_atoms=3 num_residues=1 )
626: Molecule( WAT:629 num_atoms=3 num_residues=1 )
627: Molecule( WAT:630 num_atoms=3 num_residues=1 )
628: Molecule( WAT:631 num_atoms=3 num_residues=1 )
629: Molecule( WAT:632 num_atoms=3 num_residues=1 )
)

You can also search for protein molecules using ``protein``, e.g.

>>> mols = sr.load("7SA1")
>>> print(mols["protein"])
Molecule( 7SA1:633 num_atoms=11728 num_residues=1518 )

This searches for all molecules that contain at least five residues
that are named one of the standard amino acid names.

You can see which names are matched using;

>>> print(sr.search.get_protein_residue_names())
['cyx', 'gly', 'val', 'thr', 'hip', 'his', 'tyr', 'ile', 'trp',
 'ala', 'pro', 'glh', 'ash', 'lys', 'ser', 'gln', 'arg', 'asn',
 'asp', 'cys', 'met', 'phe', 'leu', 'glu', 'hid', 'hie']

and you can see the minimum number of residues that should be
matched using

>>> print(sr.search.get_min_protein_residues())
5

You can set these parameters using the :func:`~sire.search.set_protein_residue_names`
and :func:`~sire.search.set_min_protein_residues` functions.

.. note::

   Note that the residue names are case insensitive, i.e. `ala` will
   match `ala`, `Ala`, `ALA` etc.
