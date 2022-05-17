=========
Searching
=========

Sire has a powerful search system. You have already used this to search
for atoms, residues etc by their name, number, index or element.
It is even more powerful, supporting searching by other metadata or
properties.

To explore this, let's load up the ``aladip`` system again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))

Searching by property
---------------------

You can search for atoms, residues, molecules etc by some of
their properties. Currently supported properties are mass, coordinates and charge.

Searching by mass
-----------------

For example;

>>> print(mols["atom mass < 2"])
SireMol::SelectorM<SireMol::Atom>( size=1272
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
4: MolNum(2) Atom( HA:10   [  17.29,    5.15,   16.59] )
...
1267: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1268: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1269: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1270: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1271: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

which can be shortened to

>>> print(mols["mass < 2"])
SireMol::SelectorM<SireMol::Atom>( size=1272
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
4: MolNum(2) Atom( HA:10   [  17.29,    5.15,   16.59] )
...
1267: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1268: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1269: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1270: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1271: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

will find all atoms that have a mass of less than `2 g mol-1`. You can
add the units, e.g.

>>> print(mols["mass < 2 g_per_mol"])
SireMol::SelectorM<SireMol::Atom>( size=1272
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
4: MolNum(2) Atom( HA:10   [  17.29,    5.15,   16.59] )
...
1267: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1268: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1269: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1270: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1271: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

can use any comparison you want, e.g.

>>> print(mols["mass >= 16"])
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

and also search for larger units by mass, e.g. finding all residues
that are greater than `100 g_per_mol`

>>> print(mols["residue mass > 50 g_per_mol"])
Residue( ALA:2   num_atoms=10 )

or molecules that are less than `20 g_per_mol`

>>> print(mols["molecule mass < 20 g_per_mol"])
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

or bonds where the two atoms in the bond have a total mass of greater than
25 g_per_mol

>>> print(mols["bond mass > 25 g_per_mol"])
SelectorBond( size=6
0: Bond( C:5 => O:6 )
1: Bond( C:5 => N:7 )
2: Bond( N:7 => CA:9 )
3: Bond( C:15 => O:16 )
4: Bond( C:15 => N:17 )
5: Bond( N:17 => CH3:19 )
)

Writing

>>> print(mols["mass 1.008"])
SireMol::SelectorM<SireMol::Atom>( size=1272
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
4: MolNum(2) Atom( HA:10   [  17.29,    5.15,   16.59] )
...
1267: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1268: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1269: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1270: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1271: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

is equivalent to writing

>>> print(mol["mass =~ 1.008"])
SireMol::SelectorM<SireMol::Atom>( size=1272
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
4: MolNum(2) Atom( HA:10   [  17.29,    5.15,   16.59] )
...
1267: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1268: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1269: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1270: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1271: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

where `=~` means "approximately equal to". The
`pytest algorithm <https://docs.pytest.org/en/latest/reference/reference.html#pytest-approx>`__
is used for approximate comparison. You can get the epsilon for
comparison via

>>> print(sr.search.get_approx_epsilon())
1e-06

and set it via

>>> sr.search.set_approx_epsilon(1e-6)

Searching by charge
-------------------

You can also do the same thing with charge, e.g.

>>> print(mols["charge > 0"])
SireMol::SelectorM<SireMol::Atom>( size=1275
0: MolNum(2) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(2) Atom( HH32:3  [  20.05,    3.63,   13.29] )
2: MolNum(2) Atom( HH33:4  [  18.80,    2.43,   13.73] )
3: MolNum(2) Atom( C:5     [  18.48,    4.55,   14.35] )
4: MolNum(2) Atom( H:8     [  16.68,    3.62,   14.22] )
...
1270: MolNum(630) Atom( H2:1906 [  18.40,   22.96,    9.71] )
1271: MolNum(631) Atom( H1:1908 [  22.33,    8.56,    9.83] )
1272: MolNum(631) Atom( H2:1909 [  21.07,    8.08,   10.53] )
1273: MolNum(632) Atom( H1:1911 [   8.42,   14.19,    2.16] )
1274: MolNum(632) Atom( H2:1912 [   9.90,   14.31,    1.90] )
)

gives all of the positively charged atoms, while

>>> print(mols["charge < -0.5"])
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

gives all of the atoms whose charges are less than -0.5.

The units are unit electron charges, which you can specify,

>>> print(mols["charge > 0.5 e"])
Selector<SireMol::Atom>( size=2
0:  Atom( C:5     [  18.48,    4.55,   14.35] )
1:  Atom( C:15    [  15.37,    4.19,   16.43] )
)

You can also use the same `residue`, `molecule` etc terms to search
based on the total charge on a residue, molecule etc.

>>> print(mols["residue charge 0"])
SireMol::SelectorM<SireMol::Residue>( size=631
0: MolNum(2) Residue( NME:3   num_atoms=6 )
1: MolNum(3) Residue( WAT:4   num_atoms=3 )
2: MolNum(4) Residue( WAT:5   num_atoms=3 )
3: MolNum(5) Residue( WAT:6   num_atoms=3 )
4: MolNum(6) Residue( WAT:7   num_atoms=3 )
...
626: MolNum(628) Residue( WAT:629 num_atoms=3 )
627: MolNum(629) Residue( WAT:630 num_atoms=3 )
628: MolNum(630) Residue( WAT:631 num_atoms=3 )
629: MolNum(631) Residue( WAT:632 num_atoms=3 )
630: MolNum(632) Residue( WAT:633 num_atoms=3 )
)

finds all of the neutral residues, and

>>> print(mols["bond charge < -0.5"])
SelectorBond( size=1
0: Bond( N:17 => CH3:19 )
)

finds all of the bonds where the total charge on the two atoms is
less than `-0.5 e`.

Searching by coordinates
------------------------

To search by coordinates, you can look for atoms that are within
specified distances of points or other atoms. For example,

>>> print(mols["atoms within 2.0 angstrom of element C"])
Selector<SireMol::Atom>( size=21
0:  Atom( HH31:1  [  18.45,    3.49,   12.44] )
1:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
2:  Atom( HH32:3  [  20.05,    3.63,   13.29] )
3:  Atom( HH33:4  [  18.80,    2.43,   13.73] )
4:  Atom( C:5     [  18.48,    4.55,   14.35] )
...
16:  Atom( H:18    [  15.34,    5.45,   17.96] )
17:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
18:  Atom( HH31:20 [  14.35,    3.41,   19.15] )
19:  Atom( HH32:21 [  13.19,    4.59,   18.94] )
20:  Atom( HH33:22 [  13.21,    3.33,   17.69] )
)

finds all atoms that are within `2 angstrom` of any carbon atom.
The default unit of distance is angstrom, so you could also write

>>> print(mols["atoms within 2 of element C"])
Selector<SireMol::Atom>( size=21
0:  Atom( HH31:1  [  18.45,    3.49,   12.44] )
1:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
2:  Atom( HH32:3  [  20.05,    3.63,   13.29] )
3:  Atom( HH33:4  [  18.80,    2.43,   13.73] )
4:  Atom( C:5     [  18.48,    4.55,   14.35] )
...
16:  Atom( H:18    [  15.34,    5.45,   17.96] )
17:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
18:  Atom( HH31:20 [  14.35,    3.41,   19.15] )
19:  Atom( HH32:21 [  13.19,    4.59,   18.94] )
20:  Atom( HH33:22 [  13.21,    3.33,   17.69] )
)

.. note::

    Note that we used `atoms` in this search rather than `atom`. Both are
    equivalent and can be used interchangeably. In this case, it feels better
    to write `atoms within` rather than `atom within`, but both will
    do the same thing.

You can also search by residue or other units, such as

>>> print(mols["residues within 3 angstrom of resnum 1"])
SireMol::SelectorM<SireMol::Residue>( size=14
0: MolNum(2) Residue( ACE:1   num_atoms=6 )
1: MolNum(2) Residue( ALA:2   num_atoms=10 )
2: MolNum(36) Residue( WAT:37  num_atoms=3 )
3: MolNum(55) Residue( WAT:56  num_atoms=3 )
4: MolNum(82) Residue( WAT:83  num_atoms=3 )
...
9: MolNum(466) Residue( WAT:467 num_atoms=3 )
10: MolNum(482) Residue( WAT:483 num_atoms=3 )
11: MolNum(578) Residue( WAT:579 num_atoms=3 )
12: MolNum(583) Residue( WAT:584 num_atoms=3 )
13: MolNum(623) Residue( WAT:624 num_atoms=3 )
)

returns all residues where any atom in that residue is within
`3 angstrom` of any atom in the residue with `resnum 1`.

You can also search for atoms within a point in space, e.g.

>>> print(mols["atoms within 5.0 of (0, 0, 0)"])
SireMol::SelectorM<SireMol::Atom>( size=6
0: MolNum(461) Atom( O:1397  [   1.64,    3.43,    1.98] )
1: MolNum(461) Atom( H1:1398 [   2.39,    3.63,    1.42] )
2: MolNum(461) Atom( H2:1399 [   2.00,    2.86,    2.67] )
3: MolNum(517) Atom( O:1565  [   3.69,    0.37,    0.97] )
4: MolNum(517) Atom( H1:1566 [   4.09,    1.11,    0.53] )
5: MolNum(517) Atom( H2:1567 [   3.92,    0.49,    1.90] )
)

finds all atoms within `5 angstroms` of the point `(0, 0, 0)`, while

>>> print(mols["molecules within 5.0 of (0, 0, 0)"])
SelectorMol( size=2
0: Molecule( WAT:461 num_atoms=3 num_residues=1 )
1: Molecule( WAT:517 num_atoms=3 num_residues=1 )
)

finds all molecules which have any atom that is within `5 angstroms`
of the point `(0, 0, 0)`.

Searching by molecule type
--------------------------

There are some high-level search terms that provide quick access
to searches for common types of molecules.

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

returns all water molecules. These are searched for by finding all molecules
that contain one oxygen, two hydrogens and any number of null (dummy)
atoms.

You can combine this with other searches, e.g.

>>> print(mols["water and element O"])
SireMol::SelectorM<SireMol::Atom>( size=630
0: MolNum(3) Atom( O:23    [  25.64,    8.50,   22.42] )
1: MolNum(4) Atom( O:26    [  22.83,    8.93,    4.14] )
2: MolNum(5) Atom( O:29    [   9.96,   24.84,   23.52] )
3: MolNum(6) Atom( O:32    [  25.40,   24.61,   20.94] )
4: MolNum(7) Atom( O:35    [   9.53,    4.89,   14.08] )
...
625: MolNum(628) Atom( O:1898  [  22.40,   11.84,   10.08] )
626: MolNum(629) Atom( O:1901  [   0.63,    8.69,   19.94] )
627: MolNum(630) Atom( O:1904  [  18.69,   22.12,    9.35] )
628: MolNum(631) Atom( O:1907  [  21.65,    7.88,    9.79] )
629: MolNum(632) Atom( O:1910  [   9.25,   13.73,    2.29] )
)

gives all of the oxygen atoms in water molecules.

There is a similar search term to find protein molecules.

>>> print(mols["protein"])
XXX






