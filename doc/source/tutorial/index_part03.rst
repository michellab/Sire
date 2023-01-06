=============================
Part 3 - Molecular Properties
=============================

The :class:`~sire.mol.MoleculeView`-derived classes, such as
:class:`~sire.mol.Atom`, :class:`~sire.mol.Residue` etc., are
containers for molecular information. This information is held
via properties that are associated with each view in the molecule.

To see how this works, we will first load up the ``aladip`` system.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

Properties are accessed via the ``property`` function. This takes,
as argument, the name of the property you want to retrieve.
For example, the coordinates are held in the ``coordinates`` property.

>>> print(mol.property("coordinates"))
AtomCoords( size=22
0: ( 18.4532 Å, 3.49423 Å, 12.4365 Å )
1: ( 18.9818 Å, 3.44823 Å, 13.3886 Å )
2: ( 20.0513 Å, 3.63293 Å, 13.2874 Å )
3: ( 18.798 Å, 2.43076 Å, 13.7337 Å )
4: ( 18.4805 Å, 4.54971 Å, 14.3514 Å )
...
17: ( 15.3407 Å, 5.44815 Å, 17.9626 Å )
18: ( 13.8341 Å, 3.93668 Å, 18.3509 Å )
19: ( 14.3525 Å, 3.40994 Å, 19.1521 Å )
20: ( 13.1933 Å, 4.59022 Å, 18.9428 Å )
21: ( 13.2149 Å, 3.33301 Å, 17.6874 Å )
)

and the charges are held in the ``charge`` property,

>>> print(mol.property("charge"))
SireMol::AtomCharges( size=22
0: 0.1123 |e|
1: -0.3662 |e|
2: 0.1123 |e|
3: 0.1123 |e|
4: 0.5972 |e|
...
17: 0.2719 |e|
18: -0.149 |e|
19: 0.0976 |e|
20: 0.0976 |e|
21: 0.0976 |e|
)

You can get a list of all of the properties that are contained in a molecule
via the ``property_keys`` function

>>> print(mol.property_keys())
['element', 'bond', 'dihedral', 'connectivity', 'ambertype',
 'parameters', 'gb_radius_set', 'atomtype', 'treechain', 'velocity',
 'gb_radii', 'forcefield', 'charge', 'angle', 'LJ', 'mass', 'gb_screening',
 'improper', 'intrascale', 'coordinates']

You can get all of the properties in the molecule via the
``properties`` function

>>> print(mol.properties())
Properties(
    treechain => SireMol::AtomStringProperty( size=22
0: M
1: M
2: E
3: E
4: M
...
17: E
18: M
19: E
20: E
21: E
),
    bond => TwoAtomFunctions( size=21
0:  HH31:1-CH3:2   : 340 [r - 1.09]^2
1:   CH3:2-HH32:3  : 340 [r - 1.09]^2
2:   CH3:2-HH33:4  : 340 [r - 1.09]^2
3:   CH3:2-C:5     : 317 [r - 1.522]^2
4:     C:5-O:6     : 570 [r - 1.229]^2
...
16:    N:17-H:18    : 434 [r - 1.01]^2
17:    N:17-CH3:19  : 337 [r - 1.449]^2
18:  CH3:19-HH31:20 : 340 [r - 1.09]^2
19:  CH3:19-HH32:21 : 340 [r - 1.09]^2
20:  CH3:19-HH33:22 : 340 [r - 1.09]^2
),
    dihedral => FourAtomFunctions( size=41
0:  HH31:1-CH3:2-C:5-O:6           : 0.8 cos(phi) + 0.08 cos(3 phi - 3.14159) + 0.88
1:  HH31:1-CH3:2-C:5-N:7           : 0
2:   CH3:2-C:5-N:7-H:8             : 2.5 cos(2 phi - 3.14159) + 2.5
3:   CH3:2-C:5-N:7-CA:9            : 2.5 cos(2 phi - 3.14159) + 2.5
4:  HH32:3-CH3:2-C:5-O:6           : 0.8 cos(phi) + 0.08 cos(3 phi - 3.14159) + 0.88
...
36:    O:16-C:15-N:17-H:18          : 2 cos(phi) + 2.5 cos(2 phi - 3.14159) + 4.5
37:    O:16-C:15-N:17-CH3:19        : 2.5 cos(2 phi - 3.14159) + 2.5
38:    H:18-N:17-CH3:19-HH31:20     : 0
39:    H:18-N:17-CH3:19-HH32:21     : 0
40:    H:18-N:17-CH3:19-HH33:22     : 0
),
    parameters => AmberParams( nAtoms()=22 nBonds=21, nAngles=36, nDihedrals=41 nImpropers=4 n14s=41 ),
    ambertype => SireMol::AtomStringProperty( size=22
0: HC
1: CT
2: HC
3: HC
4: C
...
17: H
18: CT
19: H1
20: H1
21: H1
),
    LJ => SireMM::AtomLJs( size=22
0: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
1: LJ( sigma = 3.39967 Å, epsilon = 0.1094 kcal mol-1 )
2: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
3: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
4: LJ( sigma = 3.39967 Å, epsilon = 0.086 kcal mol-1 )
...
17: LJ( sigma = 1.06908 Å, epsilon = 0.0157 kcal mol-1 )
18: LJ( sigma = 3.39967 Å, epsilon = 0.1094 kcal mol-1 )
19: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
20: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
21: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
),
    atomtype => SireMol::AtomStringProperty( size=22
0: HC
1: CT
2: HC
3: HC
4: C
...
17: H
18: CT
19: H1
20: H1
21: H1
),
    gb_radius_set => modified Bondi radii (mbondi)                                                   ,
    velocity => SireMol::AtomVelocities( size=22
0: ( -0.0284179 Å ps-1, -0.0279068 Å ps-1, 0.0229222 Å ps-1 )
1: ( -0.00815709 Å ps-1, 0.00301807 Å ps-1, 0.0138062 Å ps-1 )
2: ( -0.0117127 Å ps-1, 0.0275995 Å ps-1, 0.018425 Å ps-1 )
3: ( -0.0505938 Å ps-1, 0.00951136 Å ps-1, 0.0125766 Å ps-1 )
4: ( 0.0185456 Å ps-1, -0.00277097 Å ps-1, -0.005927 Å ps-1 )
...
17: ( 0.0203835 Å ps-1, -0.0236023 Å ps-1, 0.0269803 Å ps-1 )
18: ( -0.0229536 Å ps-1, 0.00122144 Å ps-1, 0.0168112 Å ps-1 )
19: ( 0.0225591 Å ps-1, -0.0584583 Å ps-1, -0.0468198 Å ps-1 )
20: ( -0.0119421 Å ps-1, 0.0177259 Å ps-1, 0.0108097 Å ps-1 )
21: ( -0.109793 Å ps-1, 0.0791381 Å ps-1, 0.0183852 Å ps-1 )
),
    charge => SireMol::AtomCharges( size=22
0: 0.1123 |e|
1: -0.3662 |e|
2: 0.1123 |e|
3: 0.1123 |e|
4: 0.5972 |e|
...
17: 0.2719 |e|
18: -0.149 |e|
19: 0.0976 |e|
20: 0.0976 |e|
21: 0.0976 |e|
),
    angle => ThreeAtomFunctions( size=36
0:  HH31:1-CH3:2-HH32:3    : 35 [theta - 1.91114]^2
1:  HH31:1-CH3:2-HH33:4    : 35 [theta - 1.91114]^2
2:  HH31:1-CH3:2-C:5       : 50 [theta - 1.91114]^2
3:   CH3:2-C:5-O:6         : 80 [theta - 2.10138]^2
4:   CH3:2-C:5-N:7         : 70 [theta - 2.03505]^2
...
31:    N:17-CH3:19-HH33:22  : 50 [theta - 1.91114]^2
32:    H:18-N:17-CH3:19     : 50 [theta - 2.06019]^2
33: HH31:20-CH3:19-HH32:21  : 35 [theta - 1.91114]^2
34: HH31:20-CH3:19-HH33:22  : 35 [theta - 1.91114]^2
35: HH32:21-CH3:19-HH33:22  : 35 [theta - 1.91114]^2
),
    coordinates => AtomCoords( size=22
0: ( 18.4532, 3.49423, 12.4365 )
1: ( 18.9818, 3.44823, 13.3886 )
2: ( 20.0513, 3.63293, 13.2874 )
3: ( 18.798, 2.43076, 13.7337 )
4: ( 18.4805, 4.54971, 14.3514 )
...
17: ( 15.3407, 5.44815, 17.9626 )
18: ( 13.8341, 3.93668, 18.3509 )
19: ( 14.3525, 3.40994, 19.1521 )
20: ( 13.1933, 4.59022, 18.9428 )
21: ( 13.2149, 3.33301, 17.6874 )
),
    element => SireMol::AtomElements( size=22
0: Hydrogen (H, 1)
1: Carbon (C, 6)
2: Hydrogen (H, 1)
3: Hydrogen (H, 1)
4: Carbon (C, 6)
...
17: Hydrogen (H, 1)
18: Carbon (C, 6)
19: Hydrogen (H, 1)
20: Hydrogen (H, 1)
21: Hydrogen (H, 1)
),
    gb_screening => SireMol::AtomFloatProperty( size=22
0: 0.85
1: 0.72
2: 0.85
3: 0.85
4: 0.72
...
17: 0.85
18: 0.72
19: 0.85
20: 0.85
21: 0.85
),
    gb_radii => SireMol::AtomRadii( size=22
0: 1.3 Å
1: 1.7 Å
2: 1.3 Å
3: 1.3 Å
4: 1.7 Å
...
17: 1.3 Å
18: 1.7 Å
19: 1.3 Å
20: 1.3 Å
21: 1.3 Å
),
    forcefield => MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine },
    intrascale => CLJNBPairs( nAtoms() == 22, nGroups() == 3 ),
    mass => SireMol::AtomMasses( size=22
0: 1.008 g mol-1
1: 12.01 g mol-1
2: 1.008 g mol-1
3: 1.008 g mol-1
4: 12.01 g mol-1
...
17: 1.008 g mol-1
18: 12.01 g mol-1
19: 1.008 g mol-1
20: 1.008 g mol-1
21: 1.008 g mol-1
),
    connectivity => Connectivity: nConnections() == 21.
Connected residues:
  * Residue ACE:1 bonded to ALA:2.
  * Residue ALA:2 bonded to ACE:1 NME:3.
  * Residue NME:3 bonded to ALA:2.
Connected atoms:
  * Atom HH31:ACE:1 bonded to CH3:ACE:1.
  * Atom CH3:ACE:1 bonded to HH31:ACE:1 C:ACE:1 HH32:ACE:1 HH33:ACE:1.
  * Atom HH32:ACE:1 bonded to CH3:ACE:1.
  * Atom HH33:ACE:1 bonded to CH3:ACE:1.
  * Atom C:ACE:1 bonded to CH3:ACE:1 N:ALA:2 O:ACE:1.
  * Atom O:ACE:1 bonded to C:ACE:1.
  * Atom N:ALA:2 bonded to CA:ALA:2 H:ALA:2 C:ACE:1.
  * Atom H:ALA:2 bonded to N:ALA:2.
  * Atom CA:ALA:2 bonded to C:ALA:2 CB:ALA:2 HA:ALA:2 N:ALA:2.
  * Atom HA:ALA:2 bonded to CA:ALA:2.
   ...,
    improper => FourAtomFunctions( size=4
0:   CH3:2-N:7-C:5-O:6             : 10.5 cos(2 phi - 3.14159) + 10.5
1:     C:5-CA:9-N:7-H:8            : 1.1 cos(2 phi - 3.14159) + 1.1
2:    CA:9-N:17-C:15-O:16          : 10.5 cos(2 phi - 3.14159) + 10.5
3:    C:15-CH3:19-N:17-H:18        : 1.1 cos(2 phi - 3.14159) + 1.1
)
)


This is a :class:`sire.base.Properties` object, which behaves a lot
like a python dictionary. For example, you can access the properties
directly from this object, e.g.

>>> print(mol.properties()["LJ"])
SireMM::AtomLJs( size=22
0: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
1: LJ( sigma = 3.39967 Å, epsilon = 0.1094 kcal mol-1 )
2: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
3: LJ( sigma = 2.64953 Å, epsilon = 0.0157 kcal mol-1 )
4: LJ( sigma = 3.39967 Å, epsilon = 0.086 kcal mol-1 )
...
17: LJ( sigma = 1.06908 Å, epsilon = 0.0157 kcal mol-1 )
18: LJ( sigma = 3.39967 Å, epsilon = 0.1094 kcal mol-1 )
19: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
20: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
21: LJ( sigma = 2.47135 Å, epsilon = 0.0157 kcal mol-1 )
)

Many of the properties shown above are atom properties. This means that
there is one value (e.g. one charge, one coordinate) per atom. You
can assign properties to any view within a molecule, e.g. residue
properties would have one value per residue, and chain properties would
have one value per chain.

Sire's property system is extremely flexible and extendable. Molecules
can have as many (or as few) properties as needed. Properties are created
from the information contained in molecular input files, but can also
be created and edited by you in a script,
e.g. :doc:`using a cursor <part03/02_cursors>`.

In this chapter you will learn how to access, use and edit properties,
across all of the views within molecules.

.. toctree::
   :maxdepth: 1

   part03/01_atom_properties
   part03/02_cursors
   part03/03_residue_properties
   part03/04_molecule_properties
   part03/05_bond_properties
   part03/06_angle_properties
   part03/07_dihedral_properties
   part03/08_improper_properties
