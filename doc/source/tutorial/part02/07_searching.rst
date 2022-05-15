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

You can search for atoms by some of their properties. Currently supported
properties are mass, coordinates and charge.

For example;

>>> print(mol["atom mass < 2"])
Selector<SireMol::Atom>( size=7440
0:  Atom( CA:2    [ -55.43,   11.35,   42.54] )
1:  Atom( C:3     [ -56.06,    9.95,   42.55] )
2:  Atom( CB:5    [ -56.32,   12.33,   41.76] )
3:  Atom( CG1:6   [ -55.68,   13.72,   41.72] )
4:  Atom( CG2:7   [ -57.70,   12.40,   42.39] )
...
7435:  Atom( C2:11705 [  10.05,   21.52,   30.16] )
7436:  Atom( C3:11706 [   9.26,   22.77,   30.56] )
7437:  Atom( C4:11708 [   8.69,   22.56,   31.97] )
7438:  Atom( C5:11709 [   7.74,   23.65,   32.46] )
7439:  Atom( C6:11712 [   8.10,   23.00,   29.56] )
)

which can be shortened to

>>> print(mol["mass < 14"])
XXX

will find all atoms that have a mass of less than `2 g mol-1`. You can
add the units, e.g.

>>> print(mol["mass < 2 g_per_mol"])
XXX

and you can use any comparison you want, e.g.

>>> print(mol["mass >= 16"])
XXX

Writing

>>> print(mol["mass 5"])
XXX

is equivalent to writing

>>> print(mol["mass == 5"])
XXX

You can also do the same thing with atomic charge, e.g.

>>> print(mol["charge > 0"])
XXX

gives all of the positively charged atoms, while

>>> print(mol["charge < -0.5"])
XXX

gives all of the atoms whose charges are less than -0.5.

The units are unit electron charges, which you can specify,

>>> print(mol["charge > 0.5 e"])
XXX

To search by coordinates, you can look for atoms that are within
specified distances of points or other atoms. For example,

>>> print(mol["atoms within 5.0 of element C"])
XXX

gives all atoms that are within 5 angstroms of any atom that matches
`element C`. You can also specify a point in space, e.g.

>>> print(mol["atoms within 3.0 of (0, 0, 0)"])
XXX

gives all atoms within 3 angstroms of the point `(0, 0, 0)`.
