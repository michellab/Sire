=============================
Dihedral Views and Properties
=============================

The ``connectivity`` property, introduced in the last chapter,
infers the presence of bonds, angles and dihedrals within a
molecule. The :class:`sire.mm.Bond`, :class:`sire.mm.Angle` and
:class:`sire.mm.Dihedral` objects provide molecule views that
can query their properties.

Dihedral views
--------------

The dihedrals in a molecule can be obtained using the ``.dihedrals()`` function,
e.g.

>>> print(mol.dihedrals())
SelectorDihedral( size=41
0: Dihedral( HH31:1 <= CH3:2 = C:5 => N:7 )
1: Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
2: Dihedral( CH3:2 <= C:5 = N:7 => CA:9 )
3: Dihedral( CH3:2 <= C:5 = N:7 => H:8 )
4: Dihedral( HH32:3 <= CH3:2 = C:5 => N:7 )
...
36: Dihedral( O:16 <= C:15 = N:17 => CH3:19 )
37: Dihedral( O:16 <= C:15 = N:17 => H:18 )
38: Dihedral( H:18 <= N:17 = CH3:19 => HH31:20 )
39: Dihedral( H:18 <= N:17 = CH3:19 => HH33:22 )
40: Dihedral( H:18 <= N:17 = CH3:19 => HH32:21 )
)

Each :class:`~sire.mm.Dihedral` object is a molecule view of the four atoms
that comprise that dihedral. It has additional functions that make it easy
to extract molecule properties that are associated with that dihedral.

For example, ``.atom0()``, ``.atom1()``, ``.atom2()`` and ``.atom3()`` can
be used to get the atoms involved in the dihedral, and from their the
coordinates of those atoms.

>>> dihedral = mol.dihedrals()[1]
>>> print(dihedral)
Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 )
>>> print(dihedral.atom0().coordinates(),
...       dihedral.atom1().coordinates(),
...       dihedral.atom2().coordinates(),
...       dihedral.atom3().coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )
( 18.4805 Å, 4.54971 Å, 14.3514 Å ) ( 19.1866 Å, 5.44143 Å, 14.7584 Å )

The dihedral object is a molecular container, so can be indexed and searched
like any other container, e.g.

>>> print(dihedral[0].coordinates(),
...       dihedral[1].coordinates(),
...       dihedral[2].coordinates(),
...       dihedral[3].coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )
( 18.4805 Å, 4.54971 Å, 14.3514 Å ) ( 19.1866 Å, 5.44143 Å, 14.7584 Å )

The ``.size()`` function is a convenience function that calculates
the size of the dihedral based on the coordinates of the atoms.

>>> print(dihedral.size())
243.281°

Note how this is reported with units. Most values calculated using sire
are returned together with their units.

You can convert to any compatible units you want using the ``.to()``
function, e.g.

>>> from sire.units import radians
>>> print(dihedral.size().to(radians))
4.246059065416177

.. note::

    The result of converting a unit is a plain double precision number.

The units system helps ensure that any calculations made in sire
make physical sense, while also reducing the risk of errors caused
by mixing units.

For example, the ``.energy()`` function on a :class:`~sire.mm.Dihedral` uses
the algebraic expressions held in the molecule's ``dihedral`` property to calculate the
energy of the dihedral. This is reported in units of kilocalories per mole.

>>> print(dihedral.energy())
0.441489 kcal mol-1

This could be converted to kilojoules per mole...

>>> from sire.units import kJ_per_mol
>>> print(dihedral.energy().to(kJ_per_mol))
1.8471895210066125

You can also get the sizes and energies of all dihedrals in a view, e.g.
to get the sizes of all dihedrals in the first residue you could use;

>>> print(mol["resnum 1"].dihedrals().sizes())
[243.281°, 5.26777°, 126.647°]

or to get the energies of all dihedrals around carbon-carbon bonds
you could use

>>> print(mol.dihedrals("element C", "element C").energies())
[0.441489 kcal mol-1, 0 , 0.217581 kcal mol-1, 0.216024 kcal mol-1,
 1.59964 kcal mol-1, 0 , 0.327286 kcal mol-1, 0 , 2.83498 kcal mol-1,
 0.756986 kcal mol-1, 0.00734048 kcal mol-1, 0.0113127 kcal mol-1,
 0.018104 kcal mol-1, 0 , 2.04634 kcal mol-1, 0 , 0 ,
 0.0619792 kcal mol-1, 0.00115906 kcal mol-1, 0.000509377 kcal mol-1,
 0.00189494 kcal mol-1, 0.00521955 kcal mol-1, 0.315842 kcal mol-1, 0 ,
 0 , 0.398546 kcal mol-1, 0.000479184 kcal mol-1, 0.0018364 kcal mol-1,
 0.00512262 kcal mol-1]

You can also use the ``.energy()`` function on a collection to get
the total energy of all dihedrals in a molecule...

>>> print(mol.dihedrals().energy())
9.53259 kcal mol-1

...or even of all dihedrals in the molecules that have been loaded
from the file.

>>> print(mols.dihedrals().energy())
9.53259 kcal mol-1

Just as for bonds, we can use a loop to find all of the dihedrals that
have a high energy, e.g.

>>> from sire.units import kcal_per_mol
>>> for dihedral in mols.dihedrals():
...     if dihedral.energy() > 0.1 * kcal_per_mol:
...         print(f"{dihedral} {dihedral.energy()}")
Dihedral( HH31:1 <= CH3:2 = C:5 => O:6 ) 0.441489 kcal mol-1
Dihedral( CH3:2 <= C:5 = N:7 => CA:9 ) 0.216024 kcal mol-1
Dihedral( CH3:2 <= C:5 = N:7 => H:8 ) 0.217581 kcal mol-1
Dihedral( HH32:3 <= CH3:2 = C:5 => O:6 ) 1.59964 kcal mol-1
Dihedral( HH33:4 <= CH3:2 = C:5 => O:6 ) 0.327286 kcal mol-1
Dihedral( C:5 <= N:7 = CA:9 => C:15 ) 0.756986 kcal mol-1
Dihedral( C:5 <= N:7 = CA:9 => CB:11 ) 2.83498 kcal mol-1
Dihedral( O:6 <= C:5 = N:7 => H:8 ) 0.115083 kcal mol-1
Dihedral( N:7 <= CA:9 = C:15 => N:17 ) 2.04634 kcal mol-1
Dihedral( HA:10 <= CA:9 = C:15 => O:16 ) 0.315842 kcal mol-1
Dihedral( CB:11 <= CA:9 = C:15 => N:17 ) 0.398546 kcal mol-1

Dihedral properties
-------------------

Just like bonds, dihedrals can also have their own per-dihedral
properties. We don't know of any molecular file formats that set
per-dihedral properties. But that doesn't stop you from setting your own!

The best way to do this is to use a cursor on the dihedral, e.g.

>>> cursor = dihedral.cursor()
>>> cursor["energy_kJ"] = dihedral.energy().to(kJ_per_mol)
>>> print(cursor["energy_kJ"])
1.84719

You can loop over lots of dihedrals to set their property, e.g.

>>> cursor = mol.cursor()
>>> for dihedral in cursor.dihedrals():
...     dihedral["energy_kJ"] = dihedral.view().energy().to(kJ_per_mol)
>>> mol = cursor.commit()
>>> print(mol.dihedrals()[1].property("energy_kJ"))
1.84719

Just for other properties, you can also use ``.apply()`` instead
of a loop.

>>> mol = mol.cursor().dihedrals().apply(
...    lambda dihedral: dihedral.set("energy_kJ", dihedral.view().energy().to(kJ_per_mol))
...   ).commit()
>>> print(mol.dihedrals()[1].property("energy_kJ"))
1.84719
