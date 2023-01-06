=========================
Bond Views and Properties
=========================

The ``connectivity`` property, introduced in the last chapter,
infers the presence of bonds, angles and dihedrals within a
molecule. The :class:`sire.mm.Bond`, :class:`sire.mm.Angle` and
:class:`sire.mm.Dihedral` objects provide molecule views that
can query their properties.

Bond views
----------

The bonds in a molecule can be obtained using the ``.bonds()`` function,
e.g.

>>> print(mol.bonds())
SelectorBond( size=21
0: Bond( HH31:1 => CH3:2 )
1: Bond( CH3:2 => HH32:3 )
2: Bond( CH3:2 => HH33:4 )
3: Bond( CH3:2 => C:5 )
4: Bond( C:5 => N:7 )
...
16: Bond( N:17 => CH3:19 )
17: Bond( N:17 => H:18 )
18: Bond( CH3:19 => HH31:20 )
19: Bond( CH3:19 => HH32:21 )
20: Bond( CH3:19 => HH33:22 )
)

Each :class:`~sire.mm.Bond` object is a molecule view of the two atoms
that comprise that bond. It has additional functions that make it easy
to extract molecule properties that are associated with that bond.

For example, ``.atom0()`` and ``.atom1()`` can be used to get the atoms
involved in the bond, and from their the coordinates of those atoms.

>>> bond = mol.bonds()[0]
>>> print(bond)
Bond( HH31:1 => CH3:2 )
>>> print(bond.atom0().coordinates(), bond.atom1().coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )

The bond object is a molecular container, so can be indexed and searched
like any other container, e.g.

>>> print(bond[0].coordinates(), bond[1].coordinates())
( 18.4532 Å, 3.49423 Å, 12.4365 Å ) ( 18.9818 Å, 3.44823 Å, 13.3886 Å )

The ``.length()`` function is a convenience function that calculates
the length of the bond based on the coordinates of the atoms.

>>> print(bond.length())
1.09 Å

Note how this is reported with units. Most values calculated using sire
are returned together with their units.

You can convert to any compatible units you want using the ``.to()``
function, e.g.

>>> from sire.units import picometer
>>> print(bond.length().to(picometer))
109.00007776856795

.. note::

    The result of converting a unit is a plain double precision number.

You can change the default units of length by calling
:func:`sire.units.set_length_unit`, e.g.

>>> sr.units.set_length_unit(picometer)
>>> print(bond.length())
109 pm

Or you can change to a full set of default "SI" units using

>>> sr.units.set_si_units()
>>> print(bond.length())
0.109 nm

You return to sire's default internal units using

>>> sr.units.set_internal_units()
>>> print(bond.length())
1.09 Å

The units system helps ensure that any calculations made in sire
make physical sense, while also reducing the risk of errors caused
by mixing units.

For example, the ``.energy()`` function on a :class:`~sire.mm.Bond` uses
the algebraic expressions held in the molecule's ``bond`` property to calculate the
energy of the bond. This is reported in units of kilocalories per mole.

>>> print(bond.energy())
2.0563e-10 kcal mol-1

This could be converted to kilojoules per mole...

>>> from sire.units import kJ_per_mol
>>> print(bond.energy().to(kJ_per_mol))
8.603571979020907e-10

but an error would be raised if you tried to convert it to picometers...

>>> print(bond.energy().to(picometer))
UserWarning: SireError::incompatible_error: Units for values
2.0563e-10 kcal mol-1 and 0.01 angstrom are incompatible.
(call sire.error.get_last_error_details() for more info)

or if you tried to add an energy to a length...

>>> print(bond.length() + bond.energy())
UserWarning: SireError::incompatible_error: Units for values 1.09 angstrom
and 2.0563e-10 kcal mol-1 are incompatible.
(call sire.error.get_last_error_details() for more info)

You can change the energy units using

>>> from sire.units import kilojoule
>>> sr.units.set_energy_unit(kilojoule)
>>> print(bond.energy())
8.60357e-10 kJ mol-1

The default energy is kilocalories, which can be reset using

>>> from sire.units import kcal
>>> sr.units.set_energy_unit(kcal)
>>> print(bond.energy())
2.0563e-10 kcal mol-1

You can also get the lengths and energies of all bonds in a view, e.g.
to get the lengths of all bonds in the first residue you could use;

>>> print(mol["resnum 1"].bonds().lengths())
[1.09 Å, 1.54643 Å, 1.09 Å, 1.09 Å, 1.20803 Å]

or to get the energies of all hydrogen-carbon bonds you
would use

>>> print(mol.bonds("element H", "element C").energies())
[2.0563e-10 kcal mol-1, 1.65144e-09 kcal mol-1, 2.2471e-13 kcal mol-1,
 7.997e-09 kcal mol-1, 1.09482e-13 kcal mol-1, 8.56699e-11 kcal mol-1,
 1.22857e-09 kcal mol-1, 2.06535e-13 kcal mol-1, 5.18497e-09 kcal mol-1,
 4.39824e-11 kcal mol-1]

You can also use the ``.energy()`` function on a collection to get
the total energy of all bonds in a molecule...

>>> print(mol.bonds().energy())
4.54821 kcal mol-1

...or even of all bonds in the molecules that have been loaded
from the file.

>>> print(mols.bonds().energy())
4.54821 kcal mol-1

This appears to be the same as the energy of the bonds in the first
molecule. We can use slicing to get the energies of all bonds except
for the first molecule.

>>> print(mols[1:].bonds().energy())
1.60207e-09 kcal mol-1

We can find the bonds that have a high energy using a loop, e.g.

>>> from sire.units import kcal_per_mol
>>> for bond in mols.bonds():
...     if bond.energy() > 0.1 * kcal_per_mol:
...         print(bond, bond.energy())
Bond( CH3:2 => C:5 ) 0.189213 kcal mol-1
Bond( C:5 => O:6 ) 0.250565 kcal mol-1
Bond( N:7 => CA:9 ) 0.27779 kcal mol-1
Bond( CA:9 => C:15 ) 0.537132 kcal mol-1
Bond( CA:9 => CB:11 ) 0.179525 kcal mol-1
Bond( C:15 => O:16 ) 0.125648 kcal mol-1
Bond( C:15 => N:17 ) 1.45641 kcal mol-1
Bond( N:17 => CH3:19 ) 1.52335 kcal mol-1

Bond properties
---------------

Just like atoms, residues and other views, bonds can also have their own per-bond
properties. Only a few molecular file formats, e.g. like the SDF format,
actually set bond properties. For example, let's load the
``cholesterol.sdf`` file.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "cholesterol.sdf"))
>>> mol = mols[0]

We can get the per-bond properties by calling the ``.property_keys()``
function.

>>> print(mol.bonds().property_keys())
['type', 'sdf_fields', 'stereoscopy']

.. note ::

    Note that the ``.bonds()`` function returns all of the bonds
    in the molecule.

You can get the value of individual properties by calling
the ``.property()`` function on a specific bond, passing in the
name of the property you want.

>>> bond = mol.bonds()[0]
>>> print(bond.property_keys())
['type', 'sdf_fields', 'stereoscopy']
>>> print(bond.property("type"))
single
>>> print(bond.property("stereoscopy"))
not stereo

.. note::

    The ``type`` property is of type :class:`sire.mol.BondType`.
    The ``stereoscopy`` property is of type :class:`sire.mol.Stereoscopy`.

You can also access the properties via a cursor on the bond, e.g.

>>> cursor = bond.cursor()
>>> print(cursor["type"])
single

You can use the cursor to edit bond properties, just like you did
for atom, residue, chain, segment and molecule properties.

>>> cursor["type"] = sr.mol.BondType.double_bond()
>>> mol = cursor.molecule().commit()
>>> print(mol.bonds()[0].property("type"))
double

You can loop over lots of bonds to set their property, e.g.

>>> cursor = mol.cursor()
>>> for bond in cursor.bonds():
...     bond["length"] = bond.view().length()
>>> mol = cursor.commit()
>>> print(mol.bonds()[1].property("length"))
1.54643 Å

Just for other properties, you can also use ``.apply()`` instead
of a loop.

>>> mol = mol.cursor().bonds().apply(
...    lambda bond: bond.set("length", bond.view().length())
...   ).commit()
>>> print(mol.bonds()[1].property("length"))
1.54643 Å
