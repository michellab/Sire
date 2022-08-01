=============================
Part 3 - Molecular Properties
=============================

The :class:`~sire.mol.MoleculeView`-derived classes, such as
:class:`~sire.mol.Atom`, :class:`~sire.mol.Residue` etc., are
containers for molecular information. This information is held
via properties that are associated with each view in the molecule.

To see how this works, we will first load up the `aladip` system.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

Properties are accessed via the `property` function. This takes,
as argument, the name of the property you want to retrieve.
For example, the coordinates are held in the `coordinates` property.

>>> print(mol.property("coordinates"))
AtomCoords( size=22
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
)

and the charges are held in the `charge` property,

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
via the `property_keys` function

>>> print(mol.property_keys())
['element', 'bond', 'dihedral', 'connectivity', 'ambertype',
 'parameters', 'gb_radius_set', 'atomtype', 'treechain', 'velocity',
 'gb_radii', 'forcefield', 'charge', 'angle', 'LJ', 'mass', 'gb_screening',
 'improper', 'intrascale', 'coordinates']

You can get all of the properties in the molecule via the
`properties` function

>>> print(mol.properties())
Properties(
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
    bond => TwoAtomFunctions( nFunctions() == 21 ),
    dihedral => FourAtomFunctions( nFunctions() == 41 ),
    connectivity => Connectivity: nConnections() == 21.
Connected residues:
  * Residue ACE:1 bonded to ALA:2.
  * Residue ALA:2 bonded to ACE:1 NME:3.
  * Residue NME:3 bonded to ALA:2.
Connected atoms:
  * Atom HH31:ACE:1 bonded to CH3:ACE:1.
  * Atom CH3:ACE:1 bonded to HH31:ACE:1 HH33:ACE:1 HH32:ACE:1 C:ACE:1.
  * Atom HH32:ACE:1 bonded to CH3:ACE:1.
  * Atom HH33:ACE:1 bonded to CH3:ACE:1.
  * Atom C:ACE:1 bonded to N:ALA:2 CH3:ACE:1 O:ACE:1.
  * Atom O:ACE:1 bonded to C:ACE:1.
  * Atom N:ALA:2 bonded to H:ALA:2 CA:ALA:2 C:ACE:1.
  * Atom H:ALA:2 bonded to N:ALA:2.
  * Atom CA:ALA:2 bonded to N:ALA:2 C:ALA:2 HA:ALA:2 CB:ALA:2.
  * Atom HA:ALA:2 bonded to CA:ALA:2.
  * Atom CB:ALA:2 bonded to HB3:ALA:2 HB2:ALA:2 CA:ALA:2 HB1:ALA:2.
  * Atom HB1:ALA:2 bonded to CB:ALA:2.
  * Atom HB2:ALA:2 bonded to CB:ALA:2.
  * Atom HB3:ALA:2 bonded to CB:ALA:2.
  * Atom C:ALA:2 bonded to N:NME:3 O:ALA:2 CA:ALA:2.
  * Atom O:ALA:2 bonded to C:ALA:2.
  * Atom N:NME:3 bonded to H:NME:3 CH3:NME:3 C:ALA:2.
  * Atom H:NME:3 bonded to N:NME:3.
  * Atom CH3:NME:3 bonded to N:NME:3 HH31:NME:3 HH33:NME:3 HH32:NME:3.
  * Atom HH31:NME:3 bonded to CH3:NME:3.
  * Atom HH32:NME:3 bonded to CH3:NME:3.
  * Atom HH33:NME:3 bonded to CH3:NME:3.,
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
    parameters => AmberParams( nAtoms()=22 nBonds=21, nAngles=36, nDihedrals=41 nImpropers=4 n14s=41 ),
    gb_radius_set => modified Bondi radii (mbondi)                                                   ,
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
    velocity => SireMol::AtomVelocities( size=22
0: ( -2.84179e-05 angstrom fs-1, -2.79068e-05 angstrom fs-1, 2.29222e-05 angstrom fs-1 )
1: ( -8.15709e-06 angstrom fs-1, 3.01807e-06 angstrom fs-1, 1.38062e-05 angstrom fs-1 )
2: ( -1.17127e-05 angstrom fs-1, 2.75995e-05 angstrom fs-1, 1.8425e-05 angstrom fs-1 )
3: ( -5.05938e-05 angstrom fs-1, 9.51136e-06 angstrom fs-1, 1.25766e-05 angstrom fs-1 )
4: ( 1.85456e-05 angstrom fs-1, -2.77097e-06 angstrom fs-1, -5.927e-06 angstrom fs-1 )
...
17: ( 2.03835e-05 angstrom fs-1, -2.36023e-05 angstrom fs-1, 2.69803e-05 angstrom fs-1 )
18: ( -2.29536e-05 angstrom fs-1, 1.22144e-06 angstrom fs-1, 1.68112e-05 angstrom fs-1 )
19: ( 2.25591e-05 angstrom fs-1, -5.84583e-05 angstrom fs-1, -4.68198e-05 angstrom fs-1 )
20: ( -1.19421e-05 angstrom fs-1, 1.77259e-05 angstrom fs-1, 1.08097e-05 angstrom fs-1 )
21: ( -0.000109793 angstrom fs-1, 7.91381e-05 angstrom fs-1, 1.83852e-05 angstrom fs-1 )
),
    gb_radii => SireMol::AtomRadii( size=22
0: 1.3 angstrom
1: 1.7 angstrom
2: 1.3 angstrom
3: 1.3 angstrom
4: 1.7 angstrom
...
17: 1.3 angstrom
18: 1.7 angstrom
19: 1.3 angstrom
20: 1.3 angstrom
21: 1.3 angstrom
),
    forcefield => MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine },
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
    angle => ThreeAtomFunctions( nFunctions() == 36 ),
    LJ => SireMM::AtomLJs( size=22
0: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
1: LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 )
2: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
3: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
4: LJ( sigma = 3.39967 A, epsilon = 0.086 kcal mol-1 )
...
17: LJ( sigma = 1.06908 A, epsilon = 0.0157 kcal mol-1 )
18: LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 )
19: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
20: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
21: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
),
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
    improper => FourAtomFunctions( nFunctions() == 4 ),
    intrascale => CLJNBPairs( nAtoms() == 22, nGroups() == 3 ),
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
)
)

This is a :class:`sire.base.Properties` object, which behaves a lot
like a python dictionary. For example, you can access the properties
directly from this object, e.g.

>>> print(mol.properties()["LJ"])
SireMM::AtomLJs( size=22
0: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
1: LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 )
2: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
3: LJ( sigma = 2.64953 A, epsilon = 0.0157 kcal mol-1 )
4: LJ( sigma = 3.39967 A, epsilon = 0.086 kcal mol-1 )
...
17: LJ( sigma = 1.06908 A, epsilon = 0.0157 kcal mol-1 )
18: LJ( sigma = 3.39967 A, epsilon = 0.1094 kcal mol-1 )
19: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
20: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
21: LJ( sigma = 2.47135 A, epsilon = 0.0157 kcal mol-1 )
)

All of the properties shown above are atom properties. This means that
there is one value (e.g. one charge, one coordinate) per atom. You
can assign properties to any view within a molecule, e.g. residue
properties would have one value per residue, and chain properties would
have one value per chain.

Sire's property system is extremely flexible and extendable. Molecules
can have as many (or as few) properties as needed. Properties are created
from the information contained in molecular input files, but can also
be created and edited by you in a script.

In this chapter you will learn how to access, use and edit properties,
across all of the views within molecules.

.. toctree::
   :maxdepth: 1

   part03/01_atom_properties
   part03/02_residue_properties
   part03/03_bond_properties
   part03/04_molecule_properties
   part04/05_editing_properties
   part04/06_cursors
