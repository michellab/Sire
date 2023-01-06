==============
Indexing Bonds
==============

Bonds represent the chemical bonds between atoms in a molecule. A
:class:`~sire.mm.Bond` is a molecular container that contains the
two atoms that make up the bond.

For example, let's look at the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

We can get all of the bonds using the :func:`~sire.mol.Molecule.bonds`
function.

>>> print(mol.bonds())
SelectorBond( size=21
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => C:5 )
2: Bond( CH3:2 => HH32:3 )
3: Bond( CH3:2 => HH33:4 )
4: Bond( C:5 => O:6 )
...
16: Bond( N:17 => H:18 )
17: Bond( N:17 => CH3:19 )
18: Bond( CH3:19 => HH32:21 )
19: Bond( CH3:19 => HH33:22 )
20: Bond( CH3:19 => HH31:20 )
)

The result (a :class:`~sire.mm.SelectorBond`) is a molecular container
for bonds. Like all molecular containers, it can be indexed,

>>> print(mol.bonds()[0])
Bond( HH31:1 => CH3:2 )

sliced

>>> print(mol.bonds()[0:5])
SelectorBond( size=5
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => C:5 )
2: Bond( CH3:2 => HH32:3 )
3: Bond( CH3:2 => HH33:4 )
4: Bond( C:5 => O:6 )
)

or accessed by a list of indicies

>>> print(mol.bonds()[ [0, 2, 4, 6, 8] ])
SelectorBond( size=5
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => C:5 )
2: Bond( C:5 => O:6 )
3: Bond( N:7 => H:8 )
4: Bond( CA:9 => HA:10 )
)

The :class:`~sire.mm.Bond` object is also a molecular container, so can
be indexed, searched and sliced just like any other container.

>>> bond = mol.bonds()[0]
>>> print(bond[0])
Atom( HH31:1  [  18.45,    3.49,   12.44] )
>>> print(bond[1])
Atom( CH3:2   [  18.98,    3.45,   13.39] )

Accessing bonds by atom
-----------------------

You can also find bonds by looking for their constituent atoms.
For example,

>>> print(mol.bonds("atomnum 1", "atomnum 2"))
SelectorBond( size=1
0: Bond( HH31:1 => CH3:2 )
)

Returns the bonds between atoms with numbers 1 and 2. If there are no
bonds that match, then an empty list is returned.

>>> print(mol.bonds("atomnum 1", "atomnum 5"))
SelectorBond::empty

If you are sure that there is only a single bond that matches, then you can use the
:func:`~sire.mol.Molecule.bond` function

>>> print(mol.bond("atomnum 1", "atomnum 2"))
Bond( HH31:1 => CH3:2 )

This will raise a ``KeyError`` if multiple bonds match, or if no bonds
match.

You can use any valid atom identifier to identify the atoms. This includes
search strings, e.g. finding all of the bonds between carbon and hydrogen
atoms.

>>> print(mol.bonds("element C", "element H"))
SelectorBond( size=10
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CA:9 => HA:10 )
4: Bond( CB:11 => HB1:12 )
5: Bond( CB:11 => HB2:13 )
6: Bond( CB:11 => HB3:14 )
7: Bond( CH3:19 => HH31:20 )
8: Bond( CH3:19 => HH32:21 )
9: Bond( CH3:19 => HH33:22 )
)

>>> print(mol.bonds("element C", "not element H"))
SelectorBond( size=9
0: Bond( CH3:2 => C:5 )
1: Bond( C:5 => O:6 )
2: Bond( C:5 => N:7 )
3: Bond( N:7 => CA:9 )
4: Bond( CA:9 => CB:11 )
5: Bond( CA:9 => C:15 )
6: Bond( C:15 => O:16 )
7: Bond( C:15 => N:17 )
8: Bond( N:17 => CH3:19 )
)

or using the atom identifying types

>>> print(mol.bonds(sr.atomid("HH31", 1), sr.atomid("CH3", 2)))
SelectorBond( size=1
0: Bond( HH31:1 => CH3:2 )
)

or using complex search strings, here finding the bonds between
atoms in two residues

>>> print(mol.bonds("atoms in residx 0", "atoms in residx 1"))
SelectorBond( size=1
0: Bond( C:5 => N:7 )
)

or mixing and matching searches

>>> print(mol.bonds(sr.atomid("C", 5), "element N"))
SelectorBond( size=1
0: Bond( C:5 => N:7 )
)

Passing in a single atom identifier will return all of the bonds
that involve that atom (or atoms).

>>> print(mol.bonds("atomnum 2"))
SelectorBond( size=4
0: Bond( CH3:2 => C:5 )
1: Bond( HH31:1 => CH3:2 )
2: Bond( CH3:2 => HH32:3 )
3: Bond( CH3:2 => HH33:4 )
)

This has returned all of the bonds that involve atom number 2, while

>>> print(mol.bonds("element C"))
SelectorBond( size=19
0: Bond( CH3:2 => HH33:4 )
1: Bond( CH3:2 => C:5 )
2: Bond( HH31:1 => CH3:2 )
3: Bond( CH3:2 => HH32:3 )
4: Bond( C:5 => O:6 )
...
14: Bond( C:15 => N:17 )
15: Bond( CH3:19 => HH31:20 )
16: Bond( CH3:19 => HH32:21 )
17: Bond( CH3:19 => HH33:22 )
18: Bond( N:17 => CH3:19 )
)

gets all of the bonds that involve carbon.

Note that you can also use ``"*"`` to match anything, so

>>> print(mol.bonds("element C", "*"))
SelectorBond( size=19
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CH3:2 => C:5 )
4: Bond( C:5 => O:6 )
...
14: Bond( C:15 => N:17 )
15: Bond( N:17 => CH3:19 )
16: Bond( CH3:19 => HH31:20 )
17: Bond( CH3:19 => HH32:21 )
18: Bond( CH3:19 => HH33:22 )
)

gives the same result.

Accessing bonds by residue
--------------------------

You can also access bonds by residue, by passing in residue identifiers.
Passing in two residues identifiers, such as here

>>> print(mol.bonds("residx 0", "residx 1"))
SelectorBond( size=1
0: Bond( C:5 => N:7 )
)

gives all of the bonds that are between those two residues.

While passing in a single residue identifier

>>> print(mol.bonds("residx 0"))
SelectorBond( size=6
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH33:4 )
2: Bond( CH3:2 => C:5 )
3: Bond( CH3:2 => HH32:3 )
4: Bond( C:5 => O:6 )
5: Bond( C:5 => N:7 )
)

gives all of the bonds that involve atoms in this residue (including the
bonds to other residues).

If you want the bonds that are contained *only* within the residue, then
use the ``bonds`` function on that residue,

>>> print(mol["residx 0"].bonds())
SelectorBond( size=5
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH33:4 )
2: Bond( CH3:2 => C:5 )
3: Bond( CH3:2 => HH32:3 )
4: Bond( C:5 => O:6 )
)

Calling the ``bonds`` function on any molecular container will return the
bonds that involve only the atoms that are fully contained in that container.

.. note::

   We have shown searching for bonds by residue. You can also search
   for bonds by chain or segment if your molecule has chains or
   segments. So ``print(mol.bonds("chainidx 0", "chainidx 1"))``
   would print the bonds between the first two chains.

Searching for bonds
-------------------

You can use search terms to look for bonds. Use ``bonds in X`` to
search for bonds within whatever matches ``X``, e.g.

>>> print(mol["bonds in residx 0"])
SelectorBond( size=5
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH33:4 )
2: Bond( CH3:2 => C:5 )
3: Bond( CH3:2 => HH32:3 )
4: Bond( C:5 => O:6 )
)

.. note::

    ``bonds in`` returns only those bonds whose atoms are wholly
    contained within whatever matches. So, in this case, these are only
    the bonds within the first residue. It doesn't include the bond
    from this residue to another residue.

If you want bonds that involve any atom in ``X`` then use
``bonds with X``, e.g.

>>> print(mol["bonds with residx 0"])
SelectorBond( size=6
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CH3:2 => C:5 )
4: Bond( C:5 => O:6 )
5: Bond( C:5 => N:7 )
)

or

>>> print(mol["bonds with atomnum 2"])
SelectorBond( size=4
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CH3:2 => C:5 )
)

or

>>> print(mol["bonds with element C"])
SelectorBond( size=19
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => C:5 )
2: Bond( CH3:2 => HH32:3 )
3: Bond( CH3:2 => HH33:4 )
4: Bond( C:5 => O:6 )
...
14: Bond( C:15 => N:17 )
15: Bond( N:17 => CH3:19 )
16: Bond( CH3:19 => HH31:20 )
17: Bond( CH3:19 => HH32:21 )
18: Bond( CH3:19 => HH33:22 )
)

You can find bonds to something using ``bonds to X``, e.g.

>>> print(mol["bonds to resnum 1"])
SelectorBond( size=1
0: Bond( C:5 => N:7 )
)

>>> print(mol["bonds to atomnum 2"])
SelectorBond( size=4
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CH3:2 => C:5 )
)

>>> print(mol["bonds to element carbon"])
SelectorBond( size=16
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( C:5 => O:6 )
4: Bond( C:5 => N:7 )
...
11: Bond( C:15 => N:17 )
12: Bond( N:17 => CH3:19 )
13: Bond( CH3:19 => HH31:20 )
14: Bond( CH3:19 => HH32:21 )
15: Bond( CH3:19 => HH33:22 )
)

.. note::

    Note that ``bonds to`` excludes bonds that are in the selection. This means
    that ``bonds to element carbon`` excludes carbon-carbon bonds. If you want
    all bonds involving carbon, then use ``bonds with element carbon``.

You can search for bonds between two groups, using
``bond from X to Y``,

>>> print(mol["bonds from resnum 1 to resnum 2"])
SelectorBond( size=1
0: Bond( C:5 => N:7 )
)

>>> print(mol["bonds from element carbon to element carbon"])
SelectorBond( size=3
0: Bond( CH3:2 => C:5 )
1: Bond( CA:9 => CB:11 )
2: Bond( CA:9 => C:15 )
)

>>> print(mol["bonds from atomnum 1 to atomnum 2"])
SelectorBond( size=1
0: Bond( HH31:1 => CH3:2 )
)

And, like all molecule containers, you can perform searches in any
container across any number of molecules, e.g.

>>> print(mols["bonds from element O to element H"])
SelectorMBond( size=1260
0: MolNum(3) Bond( O:23 => H1:24 )
1: MolNum(3) Bond( O:23 => H2:25 )
2: MolNum(4) Bond( O:26 => H1:27 )
3: MolNum(4) Bond( O:26 => H2:28 )
4: MolNum(5) Bond( O:29 => H1:30 )
...
1255: MolNum(630) Bond( O:1904 => H2:1906 )
1256: MolNum(631) Bond( O:1907 => H1:1908 )
1257: MolNum(631) Bond( O:1907 => H2:1909 )
1258: MolNum(632) Bond( O:1910 => H1:1911 )
1259: MolNum(632) Bond( O:1910 => H2:1912 )
)

.. note::

    This has returned a :class:`~sire.mm.SelectorMBond`, which is the
    multi-molecule version of the :class:`~sire.mm.SelectorBond`
    container.

Uniquely identifying a bond
---------------------------

Bonds are identified by their :class:`~sire.mol.BondID`. This is a pair
of :class:`~sire.mol.AtomID` identifiers, one for each of the two
atoms to be identified. While the atom identifier can be any type,
it is best to use atom indexes, as these uniquely identify atoms in
a molecule. A :class:`~sire.mol.BondID` comprised of two
:class:`~sire.mol.AtomIdx` identifiers will uniquely identify a single
bond.

You can easily construct a :class:`~sire.mol.BondID` using the
:func:`sire.bondid` function, e.g.

>>> print(sr.bondid(0, 1))
Bond( AtomIdx(0), AtomIdx(1) )

constructs a :class:`~sire.mol.BondID` from atom indexes,

>>> print(sr.bondid("O", "H1"))
Bond( AtomName('O'), AtomName('H1') )

constructs one from atom names, and

>>> print(sr.bondid(sr.atomid(1), sr.atomid(2)))
Bond( AtomNum(1), AtomNum(2) )

constructs one from atom numbers.

You can mix and match the IDs if you want.

You can then use the :class:`~sire.mol.BondID` to index, just like
any other identifier class.

>>> print(mols[sr.bondid("O", "H1")])
SelectorMBond( size=630
0: MolNum(3) Bond( O:23 => H1:24 )
1: MolNum(4) Bond( O:26 => H1:27 )
2: MolNum(5) Bond( O:29 => H1:30 )
3: MolNum(6) Bond( O:32 => H1:33 )
4: MolNum(7) Bond( O:35 => H1:36 )
...
625: MolNum(628) Bond( O:1898 => H1:1899 )
626: MolNum(629) Bond( O:1901 => H1:1902 )
627: MolNum(630) Bond( O:1904 => H1:1905 )
628: MolNum(631) Bond( O:1907 => H1:1908 )
629: MolNum(632) Bond( O:1910 => H1:1911 )
)

gives all of the bonds between the atoms called ``O`` and ``H1`` in
all molecules, while

>>> print(mols[0][sr.bondid(0, 1)])
Bond( HH31:1 => CH3:2 )

gives just the bond between the first and second atoms in the first molecule, and

>>> print(mols[sr.bondid(0, 1)])
SelectorMBond( size=631
0: MolNum(2) Bond( HH31:1 => CH3:2 )
1: MolNum(3) Bond( O:23 => H1:24 )
2: MolNum(4) Bond( O:26 => H1:27 )
3: MolNum(5) Bond( O:29 => H1:30 )
4: MolNum(6) Bond( O:32 => H1:33 )
...
626: MolNum(628) Bond( O:1898 => H1:1899 )
627: MolNum(629) Bond( O:1901 => H1:1902 )
628: MolNum(630) Bond( O:1904 => H1:1905 )
629: MolNum(631) Bond( O:1907 => H1:1908 )
630: MolNum(632) Bond( O:1910 => H1:1911 )
)

gives the bond between the first and second atoms in each molecule.
