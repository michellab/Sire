===============
Indexing Angles
===============

Angles represent the angle between two chemical bonds between atoms in a molecule.
There are three atoms and two bonds contained in an angle.

:class:`~sire.mm.Angle` is a molecular container that contains the
three atoms and two bonds that make up the angle.

For example, let's look at the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

We can get all of the angles using the :func:`~sire.mol.Molecule.angles`
function.

>>> print(mol.angles())
SelectorAngle( size=36
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => C:5 )
2: Angle( HH31:1 <= CH3:2 => HH33:4 )
3: Angle( CH3:2 <= C:5 => N:7 )
4: Angle( CH3:2 <= C:5 => O:6 )
...
31: Angle( N:17 <= CH3:19 => HH32:21 )
32: Angle( H:18 <= N:17 => CH3:19 )
33: Angle( HH31:20 <= CH3:19 => HH33:22 )
34: Angle( HH31:20 <= CH3:19 => HH32:21 )
35: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

The result (a :class:`~sire.mm.SelectorAngle`) is a molecular container
for angles. Like all molecular containers, it can be indexed,

>>> print(mol.angles()[0])
Angle( HH31:1 <= CH3:2 => HH32:3 )

sliced

>>> print(mol.angles()[0:5])
SelectorAngle( size=5
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => C:5 )
2: Angle( HH31:1 <= CH3:2 => HH33:4 )
3: Angle( CH3:2 <= C:5 => N:7 )
4: Angle( CH3:2 <= C:5 => O:6 )
)

or accessed by a list of indicies

>>> print(mol.angles()[ [0, 2, 4, 6, 8] ])
SelectorAngle( size=5
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( CH3:2 <= C:5 => O:6 )
3: Angle( HH32:3 <= CH3:2 => HH33:4 )
4: Angle( C:5 <= N:7 => CA:9 )
)

The :class:`~sire.mm.Angle` object is also a molecular container, so can
be indexed, searched and sliced just like any other container.

>>> angle = mol.angles()[0]
>>> print(angle[0])
Atom( HH31:1  [  18.45,    3.49,   12.44] )
>>> print(angle[1])
Atom( CH3:2   [  18.98,    3.45,   13.39] )
>>> print(angle.bonds()[0])
Bond( HH31:1 => CH3:2 )

Accessing angles by atom or bond
--------------------------------

You can also find angles by looking for their constituent atoms.
For example,

>>> print(mol.angles("atomnum 1", "atomnum 2", "atomnum 3"))
SelectorAngle( size=1
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
)

Returns the angles between atoms with numbers 1, 2 and 3. If there are no
angles that match, then an empty list is returned.

>>> print(mol.angles("atomnum 1", "atomnum 2", "atomnum 10"))
SelectorAngle::empty

If you are sure that there is only a single angle that matches, then you can use the
:func:`~sire.mol.Molecule.angle` function

>>> print(mol.angle("atomnum 1", "atomnum 2", "atomnum 3"))
Angle( HH31:1 <= CH3:2 => HH32:3 )

This will raise a ``KeyError`` if multiple angles match, or if no angles
match.

You can use any valid atom identifier to identify the atoms. This includes
search strings, e.g. here we can find all angles that have carbon
only in the center

>>> print(mol.angles("not element C", "element C", "not element C"))
SelectorAngle( size=15
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH32:3 <= CH3:2 => HH33:4 )
3: Angle( O:6 <= C:5 => N:7 )
4: Angle( N:7 <= CA:9 => HA:10 )
...
10: Angle( N:17 <= CH3:19 => HH32:21 )
11: Angle( N:17 <= CH3:19 => HH33:22 )
12: Angle( HH31:20 <= CH3:19 => HH32:21 )
13: Angle( HH31:20 <= CH3:19 => HH33:22 )
14: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

Passing in three atom identifiers, as above, will search for angles by
atom. Passing in two atom identifiers will search for angles
that contain the corresponding bond.

For example, here we can find all of the angles involving bonds
between carbon and hydrogen atoms,

>>> print(mol.angles("element C", "element H"))
SelectorAngle( size=21
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
3: Angle( HH32:3 <= CH3:2 => HH33:4 )
4: Angle( HH32:3 <= CH3:2 => C:5 )
...
16: Angle( N:17 <= CH3:19 => HH32:21 )
17: Angle( N:17 <= CH3:19 => HH33:22 )
18: Angle( HH31:20 <= CH3:19 => HH32:21 )
19: Angle( HH31:20 <= CH3:19 => HH33:22 )
20: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

This would also work using atom identifying types, e.g.
looking for angles that contains the bond between
atoms ``HH31:1`` and ``CH3:2``.

>>> print(mol.angles(sr.atomid("HH31", 1), sr.atomid("CH3", 2)))
SelectorAngle( size=3
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
)

You can even use complex search strings, here finding the angles involving
the bonds between atoms connecting two residues

>>> print(mol.angles("atoms in residx 0", "atoms in residx 1"))
SelectorAngle( size=4
0: Angle( CH3:2 <= C:5 => N:7 )
1: Angle( C:5 <= N:7 => H:8 )
2: Angle( C:5 <= N:7 => CA:9 )
3: Angle( O:6 <= C:5 => N:7 )
)

or mixing and matching searches

>>> print(mol.angles(sr.atomid("C", 5), "element N"))
SelectorAngle( size=4
0: Angle( CH3:2 <= C:5 => N:7 )
1: Angle( C:5 <= N:7 => H:8 )
2: Angle( C:5 <= N:7 => CA:9 )
3: Angle( O:6 <= C:5 => N:7 )
)

Passing in a single atom identifier will return all of the angles
that involve that atom (or atoms).

>>> print(mol.angles("atomnum 2"))
SelectorAngle( size=8
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
3: Angle( CH3:2 <= C:5 => O:6 )
4: Angle( CH3:2 <= C:5 => N:7 )
5: Angle( HH32:3 <= CH3:2 => HH33:4 )
6: Angle( HH32:3 <= CH3:2 => C:5 )
7: Angle( HH33:4 <= CH3:2 => C:5 )
)

This has returned all of the angles that involve atom number 2, while

>>> print(mol.angles("element C"))
SelectorAngle( size=36
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
3: Angle( CH3:2 <= C:5 => O:6 )
4: Angle( CH3:2 <= C:5 => N:7 )
...
31: Angle( N:17 <= CH3:19 => HH33:22 )
32: Angle( H:18 <= N:17 => CH3:19 )
33: Angle( HH31:20 <= CH3:19 => HH32:21 )
34: Angle( HH31:20 <= CH3:19 => HH33:22 )
35: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

gets all of the angles that involve carbon.

Note that you can also use ``"*"`` to match anything, so

>>> print(mol.angles("*", "element C", "*"))
SelectorAngle( size=30
0: Angle( HH31:1 <= CH3:2 => HH32:3 )
1: Angle( HH31:1 <= CH3:2 => HH33:4 )
2: Angle( HH31:1 <= CH3:2 => C:5 )
3: Angle( CH3:2 <= C:5 => O:6 )
4: Angle( CH3:2 <= C:5 => N:7 )
...
25: Angle( N:17 <= CH3:19 => HH32:21 )
26: Angle( N:17 <= CH3:19 => HH33:22 )
27: Angle( HH31:20 <= CH3:19 => HH32:21 )
28: Angle( HH31:20 <= CH3:19 => HH33:22 )
29: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

returns all of the angles that have carbon as the central atom.

Accessing angles by residue
---------------------------

You can also access angles by residue, by passing in residue identifiers.
Passing in two residue identifiers, such as here

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
