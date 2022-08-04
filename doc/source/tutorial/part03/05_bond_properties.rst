===================================
Bond, Angle and Dihedral Properties
===================================

The ``connectivity`` property, introduced in the last chapter,
infers the presence of bonds, angles and dihedrals within a
molecule. The :class:`sire.mm.Bond`, :class:`sire.mm.Angle` and
:class:`sire.mm.Dihedral` objects provide molecule views that
can query their properties.

Bond Views
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
( 18.4532, 3.49423, 12.4365 ) ( 18.9818, 3.44823, 13.3886 )

The bond object is a molecular container, so can be indexed and searched
like any other container, e.g.

>>> print(bond[0].coordinates(), bond[1].coordinates())
( 18.4532, 3.49423, 12.4365 ) ( 18.9818, 3.44823, 13.3886 )

The ``.length()`` function is a convenience function that calculates
the length of the bond based on the coordinates of the atoms.

>>> print(bond.length())
1.09 angstrom

Note how this is reported with units. Most values calculated using sire
are returned together with their units.

.. note::

    The only exception are coordinates,
    as it would be too computationally inefficient to store the units
    with the coordinates of the atoms. Coordinates are stored
    and returned in units of angstroms.

You can convert to any compatible units you want using the ``.to()``
function, e.g.

>>> from sire.units import picometer
>>> print(bond.length().to(picometer))
109.00007776856795

.. note::

    The result of converting a unit is a plain double precision number.

The units system helps ensure that any calculations made in sire
make physical sense, while also reducing the risk of errors caused
by mixing units.

For example, the ``.energy()`` function on a :class:`~sire.mm.Bond` uses
the algebraic expressions held in the ``bond`` property to calculate the
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
>>> print(mol.bonds()[0].property("length"))
1.42894 angstrom

Just for other properties, you can also use ``.apply()`` instead
of a loop.

>>> mol = mol.cursor().bonds().apply(
...    lambda bond: bond.set("length", bond.view().length())
...   ).commit()
>>> print(mol.bonds()[0].property("length"))
1.42894 angstrom


