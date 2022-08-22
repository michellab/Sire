=============================
Improper Views and Properties
=============================

The presence of bonds, angles and dihedrals are implied
by the ``connectivity`` property, and can be linked back to
real chemical features of a molecule. In contrast, impropers are
fixes to molecular mechanics forcefields that aim to either
prevent streroisomer flipping during a simulation, or to
keep the atoms around sp2-type centers (e.g. around amide bonds)
in a mostly-planar arrangement.

As such, impropers are derived from the molecular mechanics forcefield
parameters that are held in a molecule's ``improper`` property.

Improper views
--------------

The impropers defined for the forcefield of a molecule can be obtained
using the ``.impropers()`` function, e.g.

>>> print(mol.impropers())
SelectorImproper( size=4
0: Improper( CH3:2 <= C:5 => N:7 -- O:6 )
1: Improper( C:5 <= N:7 => CA:9 -- H:8 )
2: Improper( CA:9 <= C:15 => N:17 -- O:16 )
3: Improper( C:15 <= N:17 => CH3:19 -- H:18 )
)

Notice how each :class:`~sire.mm.Improper` is printed with the central atom
written as the second atom. The other three atoms are all bonded to
this central atom.

Each :class:`~sire.mm.Improper` object is a molecule view of the four atoms
that comprise that improper. It has additional functions that make it easy
to extract molecule properties that are associated with that improper.

For example, ``.atom0()``, ``.atom1()``, ``.atom2()`` and ``.atom3()`` can
be used to get the atoms involved in the improper, and from their the
coordinates of those atoms.

>>> improper = mol.impropers()[0]
>>> print(improper)
Improper( CH3:2 <= C:5 => N:7 -- O:6 )
>>> print(improper.atom0().coordinates(),
...       improper.atom1().coordinates(),
...       improper.atom2().coordinates(),
...       improper.atom3().coordinates())
( 18.9818 Å, 3.44823 Å, 13.3886 Å ) ( 18.4805 Å, 4.54971 Å, 14.3514 Å )
( 17.2214 Å, 4.31498 Å, 14.7128 Å ) ( 19.1866 Å, 5.44143 Å, 14.7584 Å )

The improper object is a molecular container, so can be indexed and searched
like any other container, e.g.

>>> print(improper[0].coordinates(),
...       improper[1].coordinates(),
...       improper[2].coordinates(),
...       improper[3].coordinates())
( 18.9818 Å, 3.44823 Å, 13.3886 Å ) ( 18.4805 Å, 4.54971 Å, 14.3514 Å )
( 17.2214 Å, 4.31498 Å, 14.7128 Å ) ( 19.1866 Å, 5.44143 Å, 14.7584 Å )

The ``.size()`` function is a convenience function that calculates
the size of the improper based on the coordinates of the atoms.

>>> print(improper.size())
-3.82426°

Note how this is reported with units. Most values calculated using sire
are returned together with their units.

You can convert to any compatible units you want using the ``.to()``
function, e.g.

>>> from sire.units import radians
>>> print(improper.size().to(radians))
-0.06674592039289437

.. note::

    The result of converting a unit is a plain double precision number.

The units system helps ensure that any calculations made in sire
make physical sense, while also reducing the risk of errors caused
by mixing units.

For example, the ``.energy()`` function on a :class:`~sire.mm.Improper` uses
the algebraic expressions held in the molecule's ``improper`` property to calculate the
energy of the improper. This is reported in units of kilocalories per mole.

>>> print(improper.energy())
0.0934184 kcal mol-1

This could be converted to kilojoules per mole...

>>> from sire.units import kJ_per_mol
>>> print(improper.energy().to(kJ_per_mol))
0.39086262536239597

You can also get the sizes and energies of all impropers in a view, e.g.
to get the sizes of all impropers in the molecule you could use;

>>> print(mol.impropers().sizes())
[-3.82426°, -0.0353552°, 4.3041°, 5.92025°]

or to get the energies of all impropers with a central nitrogen atom
you could use

>>> print(mol.impropers("*", "element N", "*", "*").energies())
[8.3952e-07 kcal mol-1, 0.0234049 kcal mol-1]

You can also use the ``.energy()`` function on a collection to get
the total energy of all impropers in a molecule...

>>> print(mol.impropers().energy())
0.235105 kcal mol-1

...or even of all impropers in the molecules that have been loaded
from the file.

>>> print(mols.impropers().energy())
0.235105 kcal mol-1

Just as for bonds, we can use a loop to find all of the impropers that
have a high energy, e.g.

>>> from sire.units import kcal_per_mol
>>> for improper in mols.impropers():
...     if improper.energy() > 0.05 * kcal_per_mol:
...         print(f"{improper} {improper.energy()}")
Improper( CH3:2 <= C:5 => N:7 -- O:6 ) 0.0934184 kcal mol-1
Improper( CA:9 <= C:15 => N:17 -- O:16 ) 0.118281 kcal mol-1

Improper properties
-------------------

Just like bonds, impropers can also have their own per-improper
properties. We don't know of any molecular file formats that set
per-improper properties. But that doesn't stop you from setting your own!

The best way to do this is to use a cursor on the improper, e.g.

>>> cursor = improper.cursor()
>>> cursor["energy_kJ"] = improper.energy().to(kJ_per_mol)
>>> print(cursor["energy_kJ"])
0.390863

You can loop over lots of impropers to set their property, e.g.

>>> cursor = mol.cursor()
>>> for improper in cursor.impropers():
...     improper["energy_kJ"] = improper.view().energy().to(kJ_per_mol)
>>> mol = cursor.commit()
>>> print(mol.impropers()[0].property("energy_kJ"))
0.390863

Just for other properties, you can also use ``.apply()`` instead
of a loop.

>>> mol = mol.cursor().impropers().apply(
...    lambda improper: improper.set("energy_kJ", improper.view().energy().to(kJ_per_mol))
...   ).commit()
>>> print(mol.impropers()[0].property("energy_kJ"))
0.390863
