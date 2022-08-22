===================
Molecule Properties
===================

Some properties apply to the molecule as a whole, and so can't be
accessed via sub-views such as atom or residue.

These properties can be accessed via the molecule. You can get a list
of properties using the `property_keys` function, e.g. loading the
`aladip` system again...

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

...we can get the list of properties in the first molecule using

>>> print(mol.property_keys())
['coordinates', 'mass', 'velocity', 'atomtype', 'bond', 'forcefield',
 'element', 'charge', 'angle', 'gb_screening', 'gb_radii', 'treechain',
 'intrascale', 'gb_radius_set', 'ambertype', 'improper', 'connectivity',
 'dihedral', 'parameters', 'LJ']

Some of these are the whole-molecule view of the sub-view properties, e.g.
accessing the ``coordinates`` property via the molecule view returns the
:class:`sire.mol.AtomCoords` object that holds all of the atomic
coordinates.

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

Other properties, such as ``connectivity`` and ``forcefield``, only make sense
from the perspective of the molecule. In this case, the `forcefield`
property holds details of the forcefield that was used to parameterise
this molecule.

>>> print(mol.property("forcefield"))
MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine }

The ``connectivity`` property holds a :class:`sire.mol.Connectivity` object
that is used to say which atoms are connected (bonded) to which
other atoms.

>>> print(mol.property("connectivity"))
Connectivity: nConnections() == 21.
Connected residues:
  * Residue ACE:1 bonded to ALA:2.
  * Residue ALA:2 bonded to ACE:1 NME:3.
  * Residue NME:3 bonded to ALA:2.
Connected atoms:
  * Atom HH31:ACE:1 bonded to CH3:ACE:1.
  * Atom CH3:ACE:1 bonded to C:ACE:1 HH31:ACE:1 HH32:ACE:1 HH33:ACE:1.
  * Atom HH32:ACE:1 bonded to CH3:ACE:1.
  * Atom HH33:ACE:1 bonded to CH3:ACE:1.
  * Atom C:ACE:1 bonded to O:ACE:1 N:ALA:2 CH3:ACE:1.
  * Atom O:ACE:1 bonded to C:ACE:1.
  * Atom N:ALA:2 bonded to C:ACE:1 H:ALA:2 CA:ALA:2.
  * Atom H:ALA:2 bonded to N:ALA:2.
  * Atom CA:ALA:2 bonded to HA:ALA:2 CB:ALA:2 N:ALA:2 C:ALA:2.
  * Atom HA:ALA:2 bonded to CA:ALA:2.
   ...

Because this molecule has been parameterised using a molecular mechanics
forcefield, it contains bond, angle and dihedral parameters.
These can be accessed via the ``bond``, ``angle`` and ``dihedral`` properties,
e.g.

>>> print(mol.property("bond"))
TwoAtomFunctions( size=21
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
)

>>> print(mol.property("angle"))
ThreeAtomFunctions( size=36
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
)

>>> print(mol.property("dihedral"))
FourAtomFunctions( size=41
0:  HH31:1-CH3:2-C:5-O:6           : 0.08 cos(3 phi - 3.14159) + 0.8 cos(phi) + 0.88
1:  HH31:1-CH3:2-C:5-N:7           : 0
2:   CH3:2-C:5-N:7-H:8             : 2.5 cos(2 phi - 3.14159) + 2.5
3:   CH3:2-C:5-N:7-CA:9            : 2.5 cos(2 phi - 3.14159) + 2.5
4:  HH32:3-CH3:2-C:5-O:6           : 0.08 cos(3 phi - 3.14159) + 0.8 cos(phi) + 0.88
...
36:    O:16-C:15-N:17-H:18          : 2.5 cos(2 phi - 3.14159) + 2 cos(phi) + 4.5
37:    O:16-C:15-N:17-CH3:19        : 2.5 cos(2 phi - 3.14159) + 2.5
38:    H:18-N:17-CH3:19-HH31:20     : 0
39:    H:18-N:17-CH3:19-HH32:21     : 0
40:    H:18-N:17-CH3:19-HH33:22     : 0
)

Instead of the forcefield parameters, the full algebraic expressions
for the bond, angle and dihedral potentials are stored. These are
stored via a `sire.cas.Expression` using sire's in-built computer
algebra system.

Complementing these, the ``intrascale`` property contains the
intramolecular non-bonded scaling factors between pairs of atoms. These
are used either to exclude atom pairs from intramolecular non-bonded
calculations, or to scale the 1-4 non-bonded interactions.

>>> print(mol.property("intrascale"))
CLJNBPairs( nAtoms() == 22, nGroups() == 3 )
