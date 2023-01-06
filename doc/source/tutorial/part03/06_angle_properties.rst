==========================
Angle Views and Properties
==========================

The ``connectivity`` property, introduced in the last chapter,
infers the presence of bonds, angles and dihedrals within a
molecule. The :class:`sire.mm.Bond`, :class:`sire.mm.Angle` and
:class:`sire.mm.Dihedral` objects provide molecule views that
can query their properties.

Angle views
-----------

The angles in a molecule can be obtained using the ``.angles()`` function,
e.g.

>>> print(mol.angles())
SelectorAngle( size=36
0: Angle( HH31:1 <= CH3:2 => C:5 )
1: Angle( HH31:1 <= CH3:2 => HH32:3 )
2: Angle( HH31:1 <= CH3:2 => HH33:4 )
3: Angle( CH3:2 <= C:5 => N:7 )
4: Angle( CH3:2 <= C:5 => O:6 )
...
31: Angle( N:17 <= CH3:19 => HH31:20 )
32: Angle( H:18 <= N:17 => CH3:19 )
33: Angle( HH31:20 <= CH3:19 => HH32:21 )
34: Angle( HH31:20 <= CH3:19 => HH33:22 )
35: Angle( HH32:21 <= CH3:19 => HH33:22 )
)

Each :class:`~sire.mm.Angle` object is a molecule view of the three atoms
that comprise that angle. It has additional functions that make it easy
to extract molecule properties that are associated with that angle.

For example, ``.atom0()``, ``.atom1()`` and ``.atom2()`` can be used to get the
atoms involved in the angle, and from their the coordinates of those atoms.

>>> angle = mol.angles()[0]
>>> print(angle)
Angle( HH31:1 <= CH3:2 => C:5 )
>>> print(angle.atom0().coordinates(),
...       angle.atom1().coordinates(),
...       angle.atom2().coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )
( 18.4805 Å, 4.54971 Å, 14.3514 Å )

The angle object is a molecular container, so can be indexed and searched
like any other container, e.g.

>>> print(angle[0].coordinates(),
...       angle[1].coordinates(),
...       angle[2].coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )
( 18.4805 Å, 4.54971 Å, 14.3514 Å )

The ``.size()`` function is a convenience function that calculates
the size of the angle based on the coordinates of the atoms.

>>> print(angle.size())
110.889°

Note how this is reported with units. Most values calculated using sire
are returned together with their units.

You can convert to any compatible units you want using the ``.to()``
function, e.g.

>>> from sire.units import radians
>>> print(angle.size().to(radians))
1.8070966851826593

.. note::

    The result of converting a unit is a plain double precision number.

The units system helps ensure that any calculations made in sire
make physical sense, while also reducing the risk of errors caused
by mixing units.

For example, the ``.energy()`` function on a :class:`~sire.mm.Angle` uses
the algebraic expressions held in the molecule's ``angle`` property to calculate the
energy of the angle. This is reported in units of kilocalories per mole.

>>> print(angle.energy())
0.378849 kcal mol-1

This could be converted to kilojoules per mole...

>>> from sire.units import kJ_per_mol
>>> print(angle.energy().to(kJ_per_mol))
1.5851034416908067

You can also get the sizes and energies of all angles in a view, e.g.
to get the sizes of all angles in the first residue you could use;

>>> print(mol["resnum 1"].angles().sizes())
[110.889°, 103.539°, 112.8°, 123.097°, 104.786°, 110.675°, 114.402°]

or to get the energies of all hydrogen-carbon-hydrogen angles you
would use

>>> print(mol.angles("element H", "element C", "element H").energies())
[0.116134 kcal mol-1, 0.378849 kcal mol-1, 0.0147076 kcal mol-1,
 0.0422992 kcal mol-1, 0.0306092 kcal mol-1, 0.125851 kcal mol-1,
 1.00296 kcal mol-1, 0.559201 kcal mol-1, 0.0010697 kcal mol-1]

You can also use the ``.energy()`` function on a collection to get
the total energy of all angles in a molecule...

>>> print(mol.angles().energy())
6.79189 kcal mol-1

...or even of all angles in the molecules that have been loaded
from the file.

>>> print(mols.angles().energy())
6.79189 kcal mol-1

Just as for bonds, we can use a loop to find all of the angles that
have a high energy, e.g.

>>> from sire.units import kcal_per_mol
>>> for angle in mols.angles():
...     if angle.energy() > 0.1 * kcal_per_mol:
...         print(f"{angle} {angle.energy()}")
Angle( HH31:1 <= CH3:2 => HH33:4 ) 0.378849 kcal mol-1
Angle( HH31:1 <= CH3:2 => HH32:3 ) 0.116134 kcal mol-1
Angle( CH3:2 <= C:5 => N:7 ) 0.794168 kcal mol-1
Angle( CH3:2 <= C:5 => O:6 ) 0.17723 kcal mol-1
Angle( HH32:3 <= CH3:2 => C:5 ) 0.338535 kcal mol-1
Angle( HH33:4 <= CH3:2 => C:5 ) 0.365931 kcal mol-1
Angle( O:6 <= C:5 => N:7 ) 0.27753 kcal mol-1
Angle( N:7 <= CA:9 => HA:10 ) 0.231417 kcal mol-1
Angle( N:7 <= CA:9 => C:15 ) 0.124129 kcal mol-1
Angle( HA:10 <= CA:9 => C:15 ) 0.105879 kcal mol-1
Angle( HB2:13 <= CB:11 => HB3:14 ) 0.125851 kcal mol-1
Angle( N:17 <= CH3:19 => HH31:20 ) 0.630416 kcal mol-1
Angle( N:17 <= CH3:19 => HH32:21 ) 0.875671 kcal mol-1
Angle( HH31:20 <= CH3:19 => HH33:22 ) 0.559201 kcal mol-1
Angle( HH31:20 <= CH3:19 => HH32:21 ) 1.00296 kcal mol-1

Angle properties
----------------

Just like bonds, angles can also have their own per-angle
properties. We don't know of any molecular file formats that set
per-angle properties. But that doesn't stop you from setting your own!

The best way to do this is to use a cursor on the angle, e.g.

>>> cursor = angle.cursor()
>>> cursor["energy_kJ"] = angle.energy().to(kJ_per_mol)
>>> print(cursor["energy_kJ"])
1.5851

You can loop over lots of angles to set their property, e.g.

>>> cursor = mol.cursor()
>>> for angle in cursor.angles():
...     angle["energy_kJ"] = angle.view().energy().to(kJ_per_mol)
>>> mol = cursor.commit()
>>> print(mol.angles()[0].property("energy_kJ"))
1.5851

Just for other properties, you can also use ``.apply()`` instead
of a loop.

>>> mol = mol.cursor().angles().apply(
...    lambda angle: angle.set("energy_kJ", angle.view().energy().to(kJ_per_mol))
...   ).commit()
>>> print(mol.angles()[0].property("energy_kJ"))
1.5851
