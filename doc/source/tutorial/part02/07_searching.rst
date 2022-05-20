=========
Searching
=========

Sire has a powerful search system. You have already used this to search
for atoms, residues etc by their name, number, index or element.
It is even more powerful, supporting searching by other metadata or
properties.

To explore this, let's load up the ``kigaki`` system.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["kigaki.gro", "kigaki.top"]))

Searching by count
------------------

You can search for all atoms in the loaded molecules using

>>> print(mols["atoms"])
SireMol::SelectorM<SireMol::Atom>( size=11120
0: MolNum(6) Atom( N:1     [  23.28,   13.14,   22.39] )
1: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
2: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
3: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
4: MolNum(6) Atom( CA:5    [  22.58,   12.65,   23.60] )
...
11115: MolNum(3622) Atom( CL:11116 [  30.22,   39.31,   33.39] )
11116: MolNum(3623) Atom( CL:11117 [  40.21,    0.35,   38.95] )
11117: MolNum(3624) Atom( CL:11118 [  42.28,   19.12,   17.40] )
11118: MolNum(3625) Atom( CL:11119 [  47.26,   45.20,   12.38] )
11119: MolNum(3626) Atom( CL:11120 [  42.51,   39.80,   21.90] )
)

or all bonds

>>> print(mols["bonds"])
SelectorBond( size=301
0: Bond( N:1 => CA:5 )
1: Bond( N:1 => H1:2 )
2: Bond( N:1 => H2:3 )
3: Bond( N:1 => H3:4 )
4: Bond( CA:5 => HA:6 )
...
296: Bond( CD:294 => HD2:296 )
297: Bond( C:298 => O:299 )
298: Bond( C:298 => N:300 )
299: Bond( N:300 => H2:302 )
300: Bond( N:300 => H1:301 )
)

(and you could do the same for `residues`, `chains`, `segments` or
`molecules`).

You can search for things that match specified numbers of atoms, residues
etc using `count( X )`. This counts the number of things within the
contained search, e.g.

>>> print(mols["count(atoms) == 3"])
SireMol::SelectorM<SireMol::Atom>( size=10797
0: MolNum(7) Atom( OW:303  [   2.30,    6.28,    1.13] )
1: MolNum(7) Atom( HW1:304 [   1.37,    6.26,    1.50] )
2: MolNum(7) Atom( HW2:305 [   2.31,    5.89,    0.21] )
3: MolNum(8) Atom( OW:306  [   2.25,    2.75,    9.96] )
4: MolNum(8) Atom( HW1:307 [   2.60,    2.58,   10.88] )
...
10792: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
10793: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
10794: MolNum(3605) Atom( OW:11097 [  37.63,   48.01,   40.24] )
10795: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
10796: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
)

This has given all of the atoms in all of the molecules where the number
of atoms was equal to 3. Typically you would want the molecules which
had this number of atoms. You can search for them using

>>> print(mols["molecules with count(atoms) == 3"])
SelectorMol( size=3599
0: Molecule( SOL:7   num_atoms=3 num_residues=1 )
1: Molecule( SOL:8   num_atoms=3 num_residues=1 )
2: Molecule( SOL:9   num_atoms=3 num_residues=1 )
3: Molecule( SOL:10  num_atoms=3 num_residues=1 )
4: Molecule( SOL:11  num_atoms=3 num_residues=1 )
...
3594: Molecule( SOL:3601 num_atoms=3 num_residues=1 )
3595: Molecule( SOL:3602 num_atoms=3 num_residues=1 )
3596: Molecule( SOL:3603 num_atoms=3 num_residues=1 )
3597: Molecule( SOL:3604 num_atoms=3 num_residues=1 )
3598: Molecule( SOL:3605 num_atoms=3 num_residues=1 )
)

Any comparison operator and any type of search is allowed, e.g. you can find all
molecules with more than 5 residues using

>>> print(mols["molecules with count(residues) > 5"])
Molecule( Protein:6 num_atoms=302 num_residues=19 )

or all residues with more than 10 atoms

THIS DOESN'T WORK, BUT IT SHOULD

>>> print(mols["residues with count(atoms) > 10"])
XXX

Searching by property
---------------------

You can search for atoms, residues, molecules etc by some of
their properties. Currently supported properties are mass, coordinates and charge.

Searching by mass
-----------------

For example;

>>> print(mols["atom mass < 2"])
SireMol::SelectorM<SireMol::Atom>( size=7370
0: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
1: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
2: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
3: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
7365: MolNum(3603) Atom( HW2:11093 [  45.51,   47.49,   46.51] )
7366: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
7367: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
7368: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
7369: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
)

which can be shortened to

>>> print(mols["mass < 2"])
SireMol::SelectorM<SireMol::Atom>( size=7370
0: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
1: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
2: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
3: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
7365: MolNum(3603) Atom( HW2:11093 [  45.51,   47.49,   46.51] )
7366: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
7367: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
7368: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
7369: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
)

will find all atoms that have a mass of less than `2 g mol-1`. You can
add the units, e.g.

>>> print(mols["mass < 2 g_per_mol"])
SireMol::SelectorM<SireMol::Atom>( size=7370
0: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
1: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
2: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
3: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
7365: MolNum(3603) Atom( HW2:11093 [  45.51,   47.49,   46.51] )
7366: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
7367: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
7368: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
7369: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
)

can use any comparison you want, e.g.

>>> print(mols["mass >= 16"])
SireMol::SelectorM<SireMol::Atom>( size=3638
0: MolNum(6) Atom( O:24    [  21.52,   14.78,   23.75] )
1: MolNum(6) Atom( O:43    [  19.51,   14.27,   26.56] )
2: MolNum(6) Atom( O:50    [  19.97,   18.72,   26.98] )
3: MolNum(6) Atom( O:60    [  22.63,   18.62,   24.67] )
4: MolNum(6) Atom( O:82    [  26.57,   19.61,   26.69] )
...
3633: MolNum(3622) Atom( CL:11116 [  30.22,   39.31,   33.39] )
3634: MolNum(3623) Atom( CL:11117 [  40.21,    0.35,   38.95] )
3635: MolNum(3624) Atom( CL:11118 [  42.28,   19.12,   17.40] )
3636: MolNum(3625) Atom( CL:11119 [  47.26,   45.20,   12.38] )
3637: MolNum(3626) Atom( CL:11120 [  42.51,   39.80,   21.90] )
)

and also search for larger units by mass, e.g. finding all residues
that are greater than `100 g_per_mol`

>>> print(mols["residue mass > 50 g_per_mol"])
Selector<SireMol::Residue>( size=18
0:  Residue( LYS:1   num_atoms=24 )
1:  Residue( ILE:2   num_atoms=19 )
2:  Residue( GLY:3   num_atoms=7 )
3:  Residue( ALA:4   num_atoms=10 )
4:  Residue( LYS:5   num_atoms=22 )
...
13:  Residue( ILE:14  num_atoms=19 )
14:  Residue( GLY:15  num_atoms=7 )
15:  Residue( ALA:16  num_atoms=10 )
16:  Residue( LYS:17  num_atoms=22 )
17:  Residue( ILE:18  num_atoms=19 )
)

or molecules that are less than `20 g_per_mol`

>>> print(mols["molecule mass < 20 g_per_mol"])
SelectorMol( size=3599
0: Molecule( SOL:7   num_atoms=3 num_residues=1 )
1: Molecule( SOL:8   num_atoms=3 num_residues=1 )
2: Molecule( SOL:9   num_atoms=3 num_residues=1 )
3: Molecule( SOL:10  num_atoms=3 num_residues=1 )
4: Molecule( SOL:11  num_atoms=3 num_residues=1 )
...
3594: Molecule( SOL:3601 num_atoms=3 num_residues=1 )
3595: Molecule( SOL:3602 num_atoms=3 num_residues=1 )
3596: Molecule( SOL:3603 num_atoms=3 num_residues=1 )
3597: Molecule( SOL:3604 num_atoms=3 num_residues=1 )
3598: Molecule( SOL:3605 num_atoms=3 num_residues=1 )
)

or bonds where the two atoms in the bond have a total mass of greater than
25 g_per_mol

>>> print(mols["bond mass > 25 g_per_mol"])
SelectorBond( size=60
0: Bond( N:1 => CA:5 )
1: Bond( CE:16 => NZ:19 )
2: Bond( C:23 => O:24 )
3: Bond( C:23 => N:25 )
4: Bond( N:25 => CA:27 )
...
55: Bond( C:279 => O:280 )
56: Bond( C:279 => N:281 )
57: Bond( N:281 => CA:283 )
58: Bond( C:298 => N:300 )
59: Bond( C:298 => O:299 )
)

Writing

>>> print(mols["mass 1.008"])
SireMol::SelectorM<SireMol::Atom>( size=7370
0: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
1: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
2: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
3: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
7365: MolNum(3603) Atom( HW2:11093 [  45.51,   47.49,   46.51] )
7366: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
7367: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
7368: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
7369: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
)

is equivalent to writing

>>> print(mol["mass =~ 1.008"])
SireMol::SelectorM<SireMol::Atom>( size=7370
0: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
1: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
2: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
3: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
7365: MolNum(3603) Atom( HW2:11093 [  45.51,   47.49,   46.51] )
7366: MolNum(3604) Atom( HW1:11095 [  41.82,   40.76,   41.82] )
7367: MolNum(3604) Atom( HW2:11096 [  41.13,   41.08,   43.27] )
7368: MolNum(3605) Atom( HW1:11098 [  38.62,   47.90,   40.15] )
7369: MolNum(3605) Atom( HW2:11099 [  37.23,   47.15,   40.56] )
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
SireMol::SelectorM<SireMol::Atom>( size=7411
0: MolNum(6) Atom( N:1     [  23.28,   13.14,   22.39] )
1: MolNum(6) Atom( H1:2    [  22.66,   13.10,   21.60] )
2: MolNum(6) Atom( H2:3    [  23.57,   14.09,   22.53] )
3: MolNum(6) Atom( H3:4    [  24.08,   12.57,   22.21] )
4: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
...
7406: MolNum(3608) Atom( NA:11102 [   9.73,   27.52,   34.34] )
7407: MolNum(3609) Atom( NA:11103 [  14.34,   30.50,   37.65] )
7408: MolNum(3610) Atom( NA:11104 [  10.83,   47.08,    0.87] )
7409: MolNum(3611) Atom( NA:11105 [  37.64,   24.06,   29.76] )
7410: MolNum(3612) Atom( NA:11106 [  45.27,   32.64,   46.48] )
)

gives all of the positively charged atoms, while

>>> print(mols["charge < -0.5"])
SireMol::SelectorM<SireMol::Atom>( size=3631
0: MolNum(6) Atom( O:24    [  21.52,   14.78,   23.75] )
1: MolNum(6) Atom( O:43    [  19.51,   14.27,   26.56] )
2: MolNum(6) Atom( O:50    [  19.97,   18.72,   26.98] )
3: MolNum(6) Atom( O:60    [  22.63,   18.62,   24.67] )
4: MolNum(6) Atom( O:82    [  26.57,   19.61,   26.69] )
...
3626: MolNum(3622) Atom( CL:11116 [  30.22,   39.31,   33.39] )
3627: MolNum(3623) Atom( CL:11117 [  40.21,    0.35,   38.95] )
3628: MolNum(3624) Atom( CL:11118 [  42.28,   19.12,   17.40] )
3629: MolNum(3625) Atom( CL:11119 [  47.26,   45.20,   12.38] )
3630: MolNum(3626) Atom( CL:11120 [  42.51,   39.80,   21.90] )
)

gives all of the atoms whose charges are less than -0.5.

The units are unit electron charges, which you can specify,

>>> print(mols["charge > 0.5 e"])
SireMol::SelectorM<SireMol::Atom>( size=25
0: MolNum(6) Atom( C:23    [  21.37,   13.56,   23.79] )
1: MolNum(6) Atom( C:42    [  19.00,   14.64,   25.50] )
2: MolNum(6) Atom( C:49    [  19.98,   17.58,   26.53] )
3: MolNum(6) Atom( C:59    [  22.96,   18.55,   25.85] )
4: MolNum(6) Atom( C:81    [  26.04,   20.38,   25.89] )
...
20: MolNum(3608) Atom( NA:11102 [   9.73,   27.52,   34.34] )
21: MolNum(3609) Atom( NA:11103 [  14.34,   30.50,   37.65] )
22: MolNum(3610) Atom( NA:11104 [  10.83,   47.08,    0.87] )
23: MolNum(3611) Atom( NA:11105 [  37.64,   24.06,   29.76] )
24: MolNum(3612) Atom( NA:11106 [  45.27,   32.64,   46.48] )
)

You can also use the same `residue`, `molecule` etc terms to search
based on the total charge on a residue, molecule etc.

>>> print(mols["residue charge 0"])
SireMol::SelectorM<SireMol::Residue>( size=3612
0: MolNum(6) Residue( ILE:2   num_atoms=19 )
1: MolNum(6) Residue( GLY:3   num_atoms=7 )
2: MolNum(6) Residue( ALA:4   num_atoms=10 )
3: MolNum(6) Residue( ILE:6   num_atoms=19 )
4: MolNum(6) Residue( ILE:8   num_atoms=19 )
...
3607: MolNum(3601) Residue( SOL:3614 num_atoms=3 )
3608: MolNum(3602) Residue( SOL:3615 num_atoms=3 )
3609: MolNum(3603) Residue( SOL:3616 num_atoms=3 )
3610: MolNum(3604) Residue( SOL:3617 num_atoms=3 )
3611: MolNum(3605) Residue( SOL:3618 num_atoms=3 )
)

finds all of the neutral residues, and

>>> print(mols["bond charge < -0.5"])
SelectorBond( size=5
0: Bond( N:61 => CA:63 )
1: Bond( N:102 => CA:104 )
2: Bond( N:160 => CA:162 )
3: Bond( N:201 => CA:203 )
4: Bond( N:259 => CA:261 )
)

finds all of the bonds where the total charge on the two atoms is
less than `-0.5 e`.

Searching by coordinates
------------------------

To search by coordinates, you can look for atoms that are within
specified distances of points or other atoms. For example,

>>> print(mols["atoms within 2.0 angstrom of element C"])
SireMol::SelectorM<SireMol::Atom>( size=268
0: MolNum(6) Atom( N:1     [  23.28,   13.14,   22.39] )
1: MolNum(6) Atom( CA:5    [  22.58,   12.65,   23.60] )
2: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
3: MolNum(6) Atom( CB:7    [  23.52,   12.72,   24.79] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
263: MolNum(6) Atom( C:298   [  20.38,   27.34,   16.01] )
264: MolNum(6) Atom( O:299   [  19.53,   26.50,   15.69] )
265: MolNum(6) Atom( N:300   [  20.91,   28.17,   15.12] )
266: MolNum(271) Atom( HW1:1096 [  17.33,   10.51,   26.59] )
267: MolNum(827) Atom( HW1:2764 [  14.66,   29.74,   19.60] )
)

finds all atoms that are within `2 angstrom` of any carbon atom.
The default unit of distance is angstrom, so you could also write

>>> print(mols["atoms within 2 of element C"])
SireMol::SelectorM<SireMol::Atom>( size=268
0: MolNum(6) Atom( N:1     [  23.28,   13.14,   22.39] )
1: MolNum(6) Atom( CA:5    [  22.58,   12.65,   23.60] )
2: MolNum(6) Atom( HA:6    [  22.30,   11.70,   23.49] )
3: MolNum(6) Atom( CB:7    [  23.52,   12.72,   24.79] )
4: MolNum(6) Atom( HB1:8   [  23.66,   13.69,   25.02] )
...
263: MolNum(6) Atom( C:298   [  20.38,   27.34,   16.01] )
264: MolNum(6) Atom( O:299   [  19.53,   26.50,   15.69] )
265: MolNum(6) Atom( N:300   [  20.91,   28.17,   15.12] )
266: MolNum(271) Atom( HW1:1096 [  17.33,   10.51,   26.59] )
267: MolNum(827) Atom( HW1:2764 [  14.66,   29.74,   19.60] )
)

.. note::

    Note that we used `atoms` in this search rather than `atom`. Both are
    equivalent and can be used interchangeably. In this case, it feels better
    to write `atoms within` rather than `atom within`, but both will
    do the same thing.

You can also search by residue or other units, such as

>>> print(mols["residues within 3 angstrom of resnum 1"])
SireMol::SelectorM<SireMol::Residue>( size=18
0: MolNum(6) Residue( LYS:1   num_atoms=24 )
1: MolNum(6) Residue( ILE:2   num_atoms=19 )
2: MolNum(6) Residue( ALA:4   num_atoms=10 )
3: MolNum(1604) Residue( SOL:1617 num_atoms=3 )
4: MolNum(1624) Residue( SOL:1637 num_atoms=3 )
...
13: MolNum(1769) Residue( SOL:1782 num_atoms=3 )
14: MolNum(1781) Residue( SOL:1794 num_atoms=3 )
15: MolNum(1800) Residue( SOL:1813 num_atoms=3 )
16: MolNum(1809) Residue( SOL:1822 num_atoms=3 )
17: MolNum(1812) Residue( SOL:1825 num_atoms=3 )
)

returns all residues where any atom in that residue is within
`3 angstrom` of any atom in the residue with `resnum 1`.

You can also search for atoms within a point in space, e.g.

>>> print(mols["atoms within 5.0 of (0, 0, 0)"])
SireMol::SelectorM<SireMol::Atom>( size=9
0: MolNum(82) Atom( OW:528  [   2.97,    0.35,    1.71] )
1: MolNum(82) Atom( HW1:529 [   3.46,    1.19,    1.50] )
2: MolNum(82) Atom( HW2:530 [   3.59,   -0.30,    2.16] )
3: MolNum(183) Atom( OW:831  [   0.75,    3.45,    0.33] )
4: MolNum(183) Atom( HW1:832 [  -0.17,    3.17,    0.04] )
5: MolNum(183) Atom( HW2:833 [   1.06,    4.22,   -0.23] )
6: MolNum(185) Atom( OW:837  [   0.72,    1.66,    3.18] )
7: MolNum(185) Atom( HW1:838 [   0.55,    2.49,    2.64] )
8: MolNum(185) Atom( HW2:839 [   1.62,    1.29,    2.96] )
)

finds all atoms within `5 angstroms` of the point `(0, 0, 0)`, while

>>> print(mols["molecules within 5.0 of (0, 0, 0)"])
SelectorMol( size=3
0: Molecule( SOL:82  num_atoms=3 num_residues=1 )
1: Molecule( SOL:183 num_atoms=3 num_residues=1 )
2: Molecule( SOL:185 num_atoms=3 num_residues=1 )
)

finds all molecules which have any atom that is within `5 angstroms`
of the point `(0, 0, 0)`.

Searching by molecule type
--------------------------

There are some high-level search terms that provide quick access
to searches for common types of molecules.

>>> print(mols["water"])
SelectorMol( size=3599
0: Molecule( SOL:7   num_atoms=3 num_residues=1 )
1: Molecule( SOL:8   num_atoms=3 num_residues=1 )
2: Molecule( SOL:9   num_atoms=3 num_residues=1 )
3: Molecule( SOL:10  num_atoms=3 num_residues=1 )
4: Molecule( SOL:11  num_atoms=3 num_residues=1 )
...
3594: Molecule( SOL:3601 num_atoms=3 num_residues=1 )
3595: Molecule( SOL:3602 num_atoms=3 num_residues=1 )
3596: Molecule( SOL:3603 num_atoms=3 num_residues=1 )
3597: Molecule( SOL:3604 num_atoms=3 num_residues=1 )
3598: Molecule( SOL:3605 num_atoms=3 num_residues=1 )
)

returns all water molecules. These are searched for by finding all molecules
that contain one oxygen, two hydrogens and any number of null (dummy)
atoms.

You can combine this with other searches, e.g.

>>> print(mols["water and element O"])
SireMol::SelectorM<SireMol::Atom>( size=3599
0: MolNum(7) Atom( OW:303  [   2.30,    6.28,    1.13] )
1: MolNum(8) Atom( OW:306  [   2.25,    2.75,    9.96] )
2: MolNum(9) Atom( OW:309  [   0.19,    3.68,    6.47] )
3: MolNum(10) Atom( OW:312  [   5.69,   12.75,   11.65] )
4: MolNum(11) Atom( OW:315  [  15.55,   15.11,    7.03] )
...
3594: MolNum(3601) Atom( OW:11085 [  49.49,   40.49,   41.73] )
3595: MolNum(3602) Atom( OW:11088 [  43.18,   44.69,   43.76] )
3596: MolNum(3603) Atom( OW:11091 [  45.83,   46.80,   45.85] )
3597: MolNum(3604) Atom( OW:11094 [  41.52,   41.48,   42.44] )
3598: MolNum(3605) Atom( OW:11097 [  37.63,   48.01,   40.24] )
)

gives all of the oxygen atoms in water molecules.

There is a similar search term to find protein molecules.

>>> print(mols["protein"])
Molecule( Protein:6 num_atoms=302 num_residues=19 )

This returns all molecules that contain at least 5 residues that have
names that are in a set of protein residue names.

You can get the set of protein residue names using;

>>> print(sr.search.get_protein_residue_names())
['hip', 'his', 'tyr', 'ile', 'trp', 'ala', 'pro', 'glh', 'ash',
 'lys', 'ser', 'gln', 'arg', 'asn', 'asp', 'cys', 'met', 'phe',
 'leu', 'glu', 'hid', 'hie', 'cyx', 'gly', 'val', 'thr']

Names are matched ignoring case, so `ALA` will be identified as a
protein residue. You can set protein residue names using

>>> sr.search.set_protein_residue_names(["ala", "ash"])
>>> print(sr.search.get_protein_residue_names())
['ala', 'ash']

You can reset the names using

>>> sr.search.set_protein_residue_names(
...     ['hip', 'his', 'tyr', 'ile', 'trp', 'ala', 'pro', 'glh', 'ash',
...      'lys', 'ser', 'gln', 'arg', 'asn', 'asp', 'cys', 'met', 'phe',
...      'leu', 'glu', 'hid', 'hie', 'cyx', 'gly', 'val', 'thr']

Similarly, you can get the minimum number of protein residues to match
using

>>> print(sr.search.get_min_protein_residues())
5

and can set it via

>>> sr.search.set_min_protein_residues(5)

Searching by custom tokens
--------------------------

It is common when searching that you will have a term that you will
want to repeat. For example, you can match proteins using `protein`,
and water molecules using `water`. You could thus match all
other molecules using

>>> print(mols["not (protein or water)"])
SelectorMol( size=21
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3607 num_atoms=1 num_residues=1 )
2: Molecule( NA:3608 num_atoms=1 num_residues=1 )
3: Molecule( NA:3609 num_atoms=1 num_residues=1 )
4: Molecule( NA:3610 num_atoms=1 num_residues=1 )
...
16: Molecule( CL:3622 num_atoms=1 num_residues=1 )
17: Molecule( CL:3623 num_atoms=1 num_residues=1 )
18: Molecule( CL:3624 num_atoms=1 num_residues=1 )
19: Molecule( CL:3625 num_atoms=1 num_residues=1 )
20: Molecule( CL:3626 num_atoms=1 num_residues=1 )
)

You can create your own search token that represents this search
using :func:`sire.search.set_token`, e.g.

>>> sr.search.set_token("other", "not (protein or water)")

This creates the token `other` that represents the search
`not (protein or water)`. You can now use `other` as a search term, e.g.

>>> print(mols["other"])
SelectorMol( size=21
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3607 num_atoms=1 num_residues=1 )
2: Molecule( NA:3608 num_atoms=1 num_residues=1 )
3: Molecule( NA:3609 num_atoms=1 num_residues=1 )
4: Molecule( NA:3610 num_atoms=1 num_residues=1 )
...
16: Molecule( CL:3622 num_atoms=1 num_residues=1 )
17: Molecule( CL:3623 num_atoms=1 num_residues=1 )
18: Molecule( CL:3624 num_atoms=1 num_residues=1 )
19: Molecule( CL:3625 num_atoms=1 num_residues=1 )
20: Molecule( CL:3626 num_atoms=1 num_residues=1 )
)

This enables you to more easily find all of the positive and negative ions,
e.g.

>>> print(mols["other and charge > 0"])
SelectorMol( size=7
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3607 num_atoms=1 num_residues=1 )
2: Molecule( NA:3608 num_atoms=1 num_residues=1 )
3: Molecule( NA:3609 num_atoms=1 num_residues=1 )
4: Molecule( NA:3610 num_atoms=1 num_residues=1 )
5: Molecule( NA:3611 num_atoms=1 num_residues=1 )
6: Molecule( NA:3612 num_atoms=1 num_residues=1 )
)

Tokens can build on one another, e.g.

>>> sr.search.set_token("positive_ions", "other and charge > 0")
>>> print(mols["positive_ions"])
SelectorMol( size=7
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3607 num_atoms=1 num_residues=1 )
2: Molecule( NA:3608 num_atoms=1 num_residues=1 )
3: Molecule( NA:3609 num_atoms=1 num_residues=1 )
4: Molecule( NA:3610 num_atoms=1 num_residues=1 )
5: Molecule( NA:3611 num_atoms=1 num_residues=1 )
6: Molecule( NA:3612 num_atoms=1 num_residues=1 )
)

You can find out what a token refers to via
:func:`sire.search.get_token`, e.g.

>>> print(sr.search.get_token("positive_ions"))
({ other => not ((protein or water)) } and charge > 0 |e|)

Note how the `other` token has been expanded into its parts.
This is because the token is expanded when it is created.
This means that the token is unaffected by what you do to
the `other` token, e.g. deleting it via

>>> sr.search.delete_token("other")

will not affect `positive_ions`

>>> print(mols["positive_ions"])
SelectorMol( size=7
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3607 num_atoms=1 num_residues=1 )
2: Molecule( NA:3608 num_atoms=1 num_residues=1 )
3: Molecule( NA:3609 num_atoms=1 num_residues=1 )
4: Molecule( NA:3610 num_atoms=1 num_residues=1 )
5: Molecule( NA:3611 num_atoms=1 num_residues=1 )
6: Molecule( NA:3612 num_atoms=1 num_residues=1 )
)

Indexing within searches
------------------------

It is often the case that multiple items will match your search. You
can request only a sub-set by indexing your search, e.g.

>>> print(mols["{positive_ions}[0]"])
Molecule( NA:3606 num_atoms=1 num_residues=1 )

returns the first item that matched the custom `positive_ions` token
you created above.

Indexing can be used with any search term. The general format is
`{search_term}[index]`. The index behaves like a python index, so can
be negative indexed or sliced, e.g.

>>> print(mols["{positive_ions}[-1]"])
Molecule( NA:3612 num_atoms=1 num_residues=1 )

>>> print(mols["{positive_ions}[0:6:2]"])
SelectorMol( size=3
0: Molecule( NA:3606 num_atoms=1 num_residues=1 )
1: Molecule( NA:3608 num_atoms=1 num_residues=1 )
2: Molecule( NA:3610 num_atoms=1 num_residues=1 )
)
