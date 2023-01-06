==================
Indexing Impropers
==================

Impropers represent the improper angle that is often used in molecular
mechanics to hold four atoms in a plane. Imagining three atoms (0, 2 and 3)
each bonded to a central atom (2), the improper angle is defined by
the four atoms 0-1-2-3 (note that there is typically no chemical bond
between atoms 2 and 3).

There are four atoms, three bonds and three angles contained in an improper.

:class:`~sire.mm.Improper` is a molecular container that contains the
four atoms, three bonds and three angles that make up the improper.

For example, let's look at the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

We can get all of the impropers using the :func:`~sire.mol.Molecule.dihedrals`
function.

>>> print(mol.impropers())
SelectorImproper( size=4
0: Improper( CA:9 <= C:15 => N:17 -- O:16 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
2: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
3: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
)

The result (a :class:`~sire.mm.SelectorImproper`) is a molecular container
for impropers. Like all molecular containers, it can be indexed,

>>> print(mol.impropers()[0])
Improper( CA:9 <= C:15 => N:17 -- O:16 )

sliced

>>> print(mol.impropers()[0:3])
SelectorImproper( size=3
0: Improper( CA:9 <= C:15 => N:17 -- O:16 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
2: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

or accessed by a list of indicies

>>> print(mol.impropers()[ [0, 2] ])
SelectorImproper( size=2
0: Improper( CA:9 <= C:15 => N:17 -- O:16 )
1: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

The :class:`~sire.mm.Improper` object is also a molecular container, so can
be indexed, searched and sliced just like any other container.

>>> improper = mol.impropers()[0]
>>> print(improper[0])
Atom( CA:9    [  16.54,    5.03,   15.81] )
>>> print(improper[1])
Atom( C:15    [  15.37,    4.19,   16.43] )
>>> print(improper.angles()[0])
Angle( CA:9 <= C:15 => N:17 )
>>> print(improper.bonds()[0])
Bond( CA:9 => C:15 )

Accessing impropers by atom, bond or angle
------------------------------------------

You can also find impropers by looking for their constituent atoms.
For example,

>>> print(mol.impropers("atomnum 9", "atomnum 15", "atomnum 17", "atomnum 16"))
SelectorImproper( size=1
0: Improper( CA:9 <= C:15 => N:17 -- O:16 )
)

Returns the impropers between atoms with numbers 9, 15, 17 and 16. If there are
no impropers that match, then an empty list is returned.

>>> print(mol.impropers("atomnum 1", "atomnum 2", "atomnum 3", "atomnum 4"))
SelectorImproper::empty

If you are sure that there is only a single improper that matches, then you can use the
:func:`~sire.mol.Molecule.improper` function

>>> print(mol.improper("atomnum 9", "atomnum 15", "atomnum 17", "atomnum 16"))
Improper( CA:9 <= C:15 => N:17 -- O:16 )

This will raise a ``KeyError`` if multiple impropers match, or if no impropers
match.

You can use any valid atom identifier to identify the atoms. This includes
search strings, e.g. here we can find all impropers nitrogen is only
found in the central (second) atom of the improper.

>>> print(mol.impropers("not element N", "element N",
...                     "not element N", "not element N"))
SelectorImproper( size=2
0: Improper( C:5 <= N:7 => CA:9 -- H:8 )
1: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

Passing in four atom identifiers, as above, will search for impropers by
atom. Passing in three atom identifiers will search for impropers that contain
the corresponding angle. For example

>>> print(mol.impropers("element C", "element N", "element C"))
SelectorImproper( size=2
0: Improper( C:5 <= N:7 => CA:9 -- H:8 )
1: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

searches for impropers that contain carbon-nitrogen-carbon angles.

Passing in two atom identifiers will search for impropers that contain
the corresponding bond.

For example, here we can find all of the impropers involving bonds
between carbon and oxygen atoms,

>>> print(mol.impropers("element C", "element O"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( CA:9 <= C:15 => N:17 -- O:16 )
)

This would also work using atom identifying types, e.g.
looking for impropers that contains the bond between
atoms ``O:6`` and ``C:5``.

>>> print(mol.impropers(sr.atomid("O", 6), sr.atomid("C", 5)))
SelectorImproper( size=1
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
)

You can even use complex search strings, here finding the impropers involving
the bonds between atoms connecting two residues

>>> print(mol.impropers("atoms in residx 0", "atoms in residx 1"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
)

or mixing and matching searches

>>> print(mol.impropers(sr.atomid("C", 5), "element N"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
)

Passing in a single atom identifier will return all of the impropers
that involve that atom (or atoms).

>>> print(mol.impropers("atomnum 5"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
)

This has returned all of the impropers that involve atom number 5, while

>>> print(mol.impropers("element H"))
SelectorImproper( size=2
0: Improper( C:5 <= N:7 => CA:9 -- H:8 )
1: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

gets all of the impropers that involve hydrogen.

Note that you can also use ``"*"`` to match anything, so

>>> print(mol.impropers("*", "element C", "*", "*"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( CA:9 <= C:15 => N:17 -- O:16 )
)

returns all of the impropers which have a carbon as the central atom.

Accessing impropers by residue
------------------------------

You can also access impropers by residue, by passing in residue identifiers.
Passing in two residue identifiers, such as here

>>> print(mol.impropers("residx 0", "residx 1"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
)

gives all of the impropers that involve bonds that are between those two residues.

While passing in a single residue identifier

>>> print(mol.impropers("residx 0"))
SelectorImproper( size=2
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
)

gives all of the impropers that involve atoms in this residue (including the
impropers to other residues).

If you want the impropers that are contained *only* within the residue, then
use the ``improperss`` function on that residue,

>>> print(mol["residx 0"].impropers())
SelectorImproper::empty

.. note::

   In this case, there are no impropers that are wholly contained
   within residues. All of the impropers in this molecule are used
   to keep the amide bonds between residues planar.

Calling the ``impropers`` function on any molecular container will return the
impropers that involve only the atoms that are fully contained in that container.

.. note::

   We have shown searching for impropers by residue. You can also search
   for impropers by chain or segment if your molecule has chains or
   segments. So ``print(mol.impropers("chainidx 0", "chainidx 1"))``
   would print the impropers between the first two chains.

Uniquely identifying an improper
--------------------------------

Impropers are identified by their :class:`~sire.mol.ImproperID`. This is a quad
of :class:`~sire.mol.AtomID` identifiers, one for each of the four
atoms to be identified. While the atom identifier can be any type,
it is best to use atom indexes, as these uniquely identify atoms in
a molecule. A :class:`~sire.mol.ImproperID` comprised of four
:class:`~sire.mol.AtomIdx` identifiers will uniquely identify a single
improper.

You can easily construct an :class:`~sire.mol.ImproperID` using the
:func:`sire.improperid` function, e.g.

>>> print(sr.improperid(0, 1, 2, 3))
Improper( AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3) )

constructs a :class:`~sire.mol.ImproperID` from atom indexes,

>>> print(sr.improperid("CH3", "C", "N", "O"))
Improper( AtomName('CH3'), AtomName('C'), AtomName('N'), AtomName('O') )

constructs one from atom names, and

>>> print(sr.improperid(sr.atomid(1), sr.atomid(2),
...                     sr.atomid(3), sr.atomid(4)))
Improper( AtomNum(1), AtomNum(2), AtomNum(3), AtomNum(4) )

constructs one from atom numbers.

You can mix and match the IDs if you want.

You can then use the :class:`~sire.mol.ImproperID` to index, just like
any other identifier class.

>>> print(mols[sr.improperid("CH3", "C", "N", "O")])
Improper( CH3:2 <= C:5 => N:7 -- O:6 )

gives the improper between the atoms called ``CH3``, ``C``, ``N`` and ``O`` in
all molecules.
