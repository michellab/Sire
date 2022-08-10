==================
Indexing Dihedrals
==================

Dihedrals represent the torsion angle around a central chemical bond between
four atoms in a molecule. There are four atoms, three bonds and two
angles contained in a dihedral.

:class:`~sire.mm.Dihedral` is a molecular container that contains the
four atoms, three bonds and two angles that make up the dihedral.

For example, let's look at the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

We can get all of the dihedrals using the :func:`~sire.mol.Molecule.dihedrals`
function.

>>> print(mol.dihedrals())
SelectorDihedral( size=41
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
2: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
3: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
4: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
...
36: Dihedral( O:16 <= C:15 = N:17 => H:18 )
37: Dihedral( O:16 <= C:15 = N:17 => CH3:19 )
38: Dihedral( H:18 <= N:17 = CH3:19 => HH32:21 )
39: Dihedral( H:18 <= N:17 = CH3:19 => HH33:22 )
40: Dihedral( H:18 <= N:17 = CH3:19 => HH31:20 )
)

The result (a :class:`~sire.mm.SelectorDihedral`) is a molecular container
for dihedrals. Like all molecular containers, it can be indexed,

>>> print(mol.dihedrals()[0])
Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )

sliced

>>> print(mol.dihedrals()[0:5])
SelectorDihedral( size=5
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
2: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
3: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
4: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
)

or accessed by a list of indicies

>>> print(mol.dihedrals()[ [0, 2, 4, 6, 8] ])
SelectorDihedral( size=5
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
2: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
3: Dihedral( HH33:4 <= CH3:2 = C:5 => N:7 )
4: Dihedral( C:5 <= N:7 = CA:9 => C:15 )
)

The :class:`~sire.mm.Dihedral` object is also a molecular container, so can
be indexed, searched and sliced just like any other container.

>>> dihedral = mol.dihedrals()[0]
>>> print(dihedral[0])
Atom( HH31:1  [  18.45,    3.49,   12.44] )
>>> print(dihedral[1])
Atom( CH3:2   [  18.98,    3.45,   13.39] )
>>> print(dihedral.angles()[0])
Angle( HH31:1 <= CH3:2 => C:5 )
>>> print(dihedral.bonds()[0])
Bond( HH31:1 => CH3:2 )

Accessing dihedrals by atom, bond or angle
------------------------------------------

You can also find dihedrals by looking for their constituent atoms.
For example,

>>> print(mol.dihedrals("atomnum 1", "atomnum 2", "atomnum 5", "atomnum 7"))
SelectorDihedral( size=1
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
)

Returns the dihedrals between atoms with numbers 1, 2, 5 and 7. If there are no
dihedrals that match, then an empty list is returned.

>>> print(mol.dihedrals("atomnum 1", "atomnum 2", "atomnum 3", "atomnum 4"))
SelectorDihedral::empty

If you are sure that there is only a single dihedral that matches, then you can use the
:func:`~sire.mol.Molecule.dihedral` function

>>> print(mol.dihedral("atomnum 1", "atomnum 2", "atomnum 5", "atomnum 7"))
Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )

This will raise a ``KeyError`` if multiple dihedrals match, or if no dihedrals
match.

You can use any valid atom identifier to identify the atoms. This includes
search strings, e.g. here we can find all dihedrals where the carbon atoms
are only in the middle two atoms.

>>> print(mol.dihedrals("not element C", "element C",
...                     "element C", "not element C"))
SelectorDihedral( size=16
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
2: Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 )
3: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
4: Dihedral( HH33:4 <= CH3:2 = C:5 => O:6 )
...
11: Dihedral( HA:10 <= CA:9 = CB:11 => HB1:12 )
12: Dihedral( HA:10 <= CA:9 = CB:11 => HB2:13 )
13: Dihedral( HA:10 <= CA:9 = CB:11 => HB3:14 )
14: Dihedral( HA:10 <= CA:9 = C:15 => O:16 )
15: Dihedral( HA:10 <= CA:9 = C:15 => N:17 )
)

Passing in four atom identifiers, as above, will search for dihedrals by
atom. Passing in three atom identifiers will search for dihedrals that contain
the corresponding angle. For example

>>> print(mol.dihedrals("element H", "element C", "element C"))
SelectorDihedral( size=17
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
2: Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 )
3: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
4: Dihedral( HH33:4 <= CH3:2 = C:5 => O:6 )
...
12: Dihedral( HA:10 <= CA:9 = C:15 => O:16 )
13: Dihedral( HA:10 <= CA:9 = C:15 => N:17 )
14: Dihedral( HB1:12 <= CB:11 = CA:9 => C:15 )
15: Dihedral( HB2:13 <= CB:11 = CA:9 => C:15 )
16: Dihedral( HB3:14 <= CB:11 = CA:9 => C:15 )
)

searches for dihedrals that contain angles between hydrogen-carbon-carbon
atoms.

Passing in two atom identifiers will search for dihedrals that contain
the corresponding bond.

For example, here we can find all of the dihedrals involving bonds
between carbon and hydrogen atoms,

>>> print(mol.dihedrals("element C", "element H"))
SelectorDihedral( size=25
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
2: Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 )
3: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
4: Dihedral( HH33:4 <= CH3:2 = C:5 => O:6 )
...
20: Dihedral( C:15 <= N:17 = CH3:19 => HH32:21 )
21: Dihedral( C:15 <= N:17 = CH3:19 => HH33:22 )
22: Dihedral( H:18 <= N:17 = CH3:19 => HH31:20 )
23: Dihedral( H:18 <= N:17 = CH3:19 => HH32:21 )
24: Dihedral( H:18 <= N:17 = CH3:19 => HH33:22 )
)

This would also work using atom identifying types, e.g.
looking for dihedrals that contains the bond between
atoms ``HH31:1`` and ``CH3:2``.

>>> print(mol.dihedrals(sr.atomid("HH31", 1), sr.atomid("CH3", 2)))
SelectorDihedral( size=2
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
)

You can even use complex search strings, here finding the dihedrals involving
the bonds between atoms connecting two residues

>>> print(mol.dihedrals("atoms in residx 0", "atoms in residx 1"))
SelectorDihedral( size=10
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
2: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
3: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
4: Dihedral( HH33:4 <= CH3:2 = C:5 => N:7 )
5: Dihedral( C:5 <= N:7 = CA:9 => HA:10 )
6: Dihedral( C:5 <= N:7 = CA:9 => CB:11 )
7: Dihedral( C:5 <= N:7 = CA:9 => C:15 )
8: Dihedral( O:6 <= C:5 = N:7 => H:8 )
9: Dihedral( O:6 <= C:5 = N:7 => CA:9 )
)

or mixing and matching searches

>>> print(mol.dihedrals(sr.atomid("C", 5), "element N"))
SelectorDihedral( size=10
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
2: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
3: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
4: Dihedral( HH33:4 <= CH3:2 = C:5 => N:7 )
5: Dihedral( C:5 <= N:7 = CA:9 => HA:10 )
6: Dihedral( C:5 <= N:7 = CA:9 => CB:11 )
7: Dihedral( C:5 <= N:7 = CA:9 => C:15 )
8: Dihedral( O:6 <= C:5 = N:7 => H:8 )
9: Dihedral( O:6 <= C:5 = N:7 => CA:9 )
)

Passing in a single atom identifier will return all of the dihedrals
that involve that atom (or atoms).

>>> print(mol.dihedrals("atomnum 2"))
SelectorDihedral( size=8
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
2: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
3: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
4: Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 )
5: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
6: Dihedral( HH33:4 <= CH3:2 = C:5 => O:6 )
7: Dihedral( HH33:4 <= CH3:2 = C:5 => N:7 )
)

This has returned all of the dihedrals that involve atom number 2, while

>>> print(mol.dihedrals("element C"))
SelectorDihedral( size=41
0: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
2: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
3: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
4: Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 )
...
36: Dihedral( O:16 <= C:15 = N:17 => H:18 )
37: Dihedral( O:16 <= C:15 = N:17 => CH3:19 )
38: Dihedral( H:18 <= N:17 = CH3:19 => HH31:20 )
39: Dihedral( H:18 <= N:17 = CH3:19 => HH32:21 )
40: Dihedral( H:18 <= N:17 = CH3:19 => HH33:22 )
)

gets all of the dihedrals that involve carbon.

Note that you can also use ``"*"`` to match anything, so

>>> print(mol.dihedrals("*", "element C", "element N", "*"))
SelectorDihedral( size=20
0: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
1: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
2: Dihedral( C:5 <= N:7 = CA:9 => HA:10 )
3: Dihedral( C:5 <= N:7 = CA:9 => CB:11 )
4: Dihedral( C:5 <= N:7 = CA:9 => C:15 )
...
15: Dihedral( O:16 <= C:15 = N:17 => H:18 )
16: Dihedral( O:16 <= C:15 = N:17 => CH3:19 )
17: Dihedral( H:18 <= N:17 = CH3:19 => HH31:20 )
18: Dihedral( H:18 <= N:17 = CH3:19 => HH32:21 )
19: Dihedral( H:18 <= N:17 = CH3:19 => HH33:22 )
)

returns all of the dihedrals that are around carbon-nitrogen bonds.

Accessing angles by residue
---------------------------

You can also access angles by residue, by passing in residue identifiers.
Passing in two residues identifiers, such as here

>>> print(mol.angles("residx 0", "residx 1"))
SelectorAngle( size=4
0: Angle( CH3:2 <= C:5 => N:7 )
1: Angle( C:5 <= N:7 => H:8 )
2: Angle( C:5 <= N:7 => CA:9 )
3: Angle( O:6 <= C:5 => N:7 )
)

gives all of the angles that involve bonds that are between those two residues.

While passing in a single residue identifier

>>> print(mol.angles("residx 0"))
SelectorAngle( size=11
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
3: Angle( CH3:2 <= C:5 => O:6 )
4: Angle( CH3:2 <= C:5 => N:7 )
...
6: Angle( HH32:3 <= CH3:2 => C:5 )
7: Angle( HH33:4 <= CH3:2 => C:5 )
8: Angle( C:5 <= N:7 => H:8 )
9: Angle( C:5 <= N:7 => CA:9 )
10: Angle( O:6 <= C:5 => N:7 )
)

gives all of the angles that involve atoms in this residue (including the
angles to other residues).

If you want the angles that are contained *only* within the residue, then
use the ``angles`` function on that residue,

>>> print(mol["residx 0"].angles())
SelectorAngle( size=7
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => C:5 )
2: Angle( HH31:1 <= CH3:2 => HH33:4 )
3: Angle( CH3:2 <= C:5 => O:6 )
4: Angle( HH32:3 <= CH3:2 => C:5 )
5: Angle( HH32:3 <= CH3:2 => HH33:4 )
6: Angle( HH33:4 <= CH3:2 => C:5 )
)

Calling the ``angles`` function on any molecular container will return the
angles that involve only the atoms that are fully contained in that container.

.. note::

   We have shown searching for angles by residue. You can also search
   for angles by chain or segment if your molecule has chains or
   segments. So ``print(mol.angles("chainidx 0", "chainidx 1"))``
   would print the angles between the first two chains.

Uniquely identifying an angle
-----------------------------

Angles are identified by their :class:`~sire.mol.AngleID`. This is a triple
of :class:`~sire.mol.AtomID` identifiers, one for each of the three
atoms to be identified. While the atom identifier can be any type,
it is best to use atom indexes, as these uniquely identify atoms in
a molecule. A :class:`~sire.mol.AngleID` comprised of three
:class:`~sire.mol.AtomIdx` identifiers will uniquely identify a single
angle.

You can easily construct a :class:`~sire.mol.AngleID` using the
:func:`sire.angleid` function, e.g.

>>> print(sr.angleid(0, 1, 2))
Angle( AtomIdx(0), AtomIdx(1), AtomIdx(2) )

constructs a :class:`~sire.mol.AngleID` from atom indexes,

>>> print(sr.angleid("H2", "O", "H1"))
Angle( AtomName('H2'), AtomName('O'), AtomName('H1') )

constructs one from atom names, and

>>> print(sr.angleid(sr.atomid(1), sr.atomid(2), sr.atomid(3)))
Angle( AtomNum(1), AtomNum(2), AtomNum(3) )

constructs one from atom numbers.

You can mix and match the IDs if you want.

You can then use the :class:`~sire.mol.AngleID` to index, just like
any other identifier class.

>>> print(mols[sr.angleid("H2", "O", "H1")])
SelectorMAngle( size=630
0: MolNum(3) Angle( H1:24 <= O:23 => H2:25 )
1: MolNum(4) Angle( H1:27 <= O:26 => H2:28 )
2: MolNum(5) Angle( H1:30 <= O:29 => H2:31 )
3: MolNum(6) Angle( H1:33 <= O:32 => H2:34 )
4: MolNum(7) Angle( H1:36 <= O:35 => H2:37 )
...
625: MolNum(628) Angle( H1:1899 <= O:1898 => H2:1900 )
626: MolNum(629) Angle( H1:1902 <= O:1901 => H2:1903 )
627: MolNum(630) Angle( H1:1905 <= O:1904 => H2:1906 )
628: MolNum(631) Angle( H1:1908 <= O:1907 => H2:1909 )
629: MolNum(632) Angle( H1:1911 <= O:1910 => H2:1912 )
)

gives all of the angles between the atoms called ``H2``, ``O`` and ``H1`` in
all molecules, while

>>> print(mols[0][sr.angleid(0, 1, 2)])
Angle( HH31:1 <= CH3:2 => HH32:3 )

gives just the angle between the first three atoms in the first molecule, and

>>> print(mols[sr.angleid(0, 1, 2)])
SelectorMAngle( size=631
0: MolNum(2) Angle( HH31:1 <= CH3:2 => HH32:3 )
1: MolNum(3) Angle( O:23 <= H1:24 => H2:25 )
2: MolNum(4) Angle( O:26 <= H1:27 => H2:28 )
3: MolNum(5) Angle( O:29 <= H1:30 => H2:31 )
4: MolNum(6) Angle( O:32 <= H1:33 => H2:34 )
...
626: MolNum(628) Angle( O:1898 <= H1:1899 => H2:1900 )
627: MolNum(629) Angle( O:1901 <= H1:1902 => H2:1903 )
628: MolNum(630) Angle( O:1904 <= H1:1905 => H2:1906 )
629: MolNum(631) Angle( O:1907 <= H1:1908 => H2:1909 )
630: MolNum(632) Angle( O:1910 <= H1:1911 => H2:1912 )
)

gives the angle between the first three atoms in each molecule.
