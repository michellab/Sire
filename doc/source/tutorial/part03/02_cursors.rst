==============================
Cursors and Editing Properties
==============================

It is easy to edit the properties of a molecule. The best way to
do this is by creating a :class:`sire.mol.Cursor` for the molecule.

>>> cursor = mol.cursor()
>>> print(cursor)
Cursor(molecule, ACE:2)

The :class:`~sire.mol.Cursor` represents a property editor that can
be scanned across the molecule. First, we will scan to the atom
called ``CA`` using the ``atom`` function.

>>> cursor = cursor.atom("CA")
>>> print(cursor)
Cursor(atom, CA:9)

The cursor provides a dictionary-style interface to all of the atom's
properties.

>>> print(cursor.keys())
['gb_screening', 'mass', 'ambertype', 'coordinates',
 'treechain', 'atomtype', 'velocity', 'LJ', 'gb_radii', 'charge', 'element']
>>> print(cursor.items())
[('gb_radii', 1.7 Å),
 ('coordinates', ( 16.5371 Å, 5.02707 Å, 15.812 Å )),
 ('charge', 0.0337 |e|),
 ('mass', 12.01 g mol-1),
 ('LJ', LJ( sigma = 3.39967 Å, epsilon = 0.1094 kcal mol-1 )),
 ('velocity', ( 0.00308134 Å ps-1, -0.0190426 Å ps-1, 0.00618047 Å ps-1 )),
 ('atomtype', 'CX'),
 ('gb_screening', 0.72),
 ('treechain', 'M   '),
 ('element', Carbon (C, 6)),
 ('ambertype', 'CX')]
>>> print(cursor["coordinates"])
( 16.5371 Å, 5.02707 Å, 15.812 Å )

Assiging to the dictionary will update the corresponding property.

>>> cursor["coordinates"] = (1.0, 2.0, 3.0)
>>> print(cursor["coordinates"])
( 1 Å, 2 Å, 3 Å )

Assigning to a non-existent key will create a new property.

>>> cursor["color"] = "charcoal"
>>> print(cursor["color"])
charcoal

An alternative to using the dictionary-type functions is to use
the ``set`` and ``get`` functions, e.g.

>>> cursor.set("color", "charcoal")
>>> print(cursor.get("color"))
charcoal

The cursor is editing a copy of the molecule. To commit and save the
changes, you need to use the ``commit()`` function.

>>> mol = cursor.molecule().commit()
>>> print(mol["CA"].property("color"))
charcoal
>>> print(mol["CA"].coordinates())
( 1 Å, 2 Å, 3 Å )

.. note::

    The ``mol`` object was itself a copy from the original in the
    ``mols`` container loaded from the file. To update the original
    in ``mols``, you need to call the ``update`` function, e.g.
    ``mols.update(mol)``. We will cover editing and updating
    molecules in more detail in a later chapter.

Cursors are useful as they make it easy to iterate over and edit
the properties of several atoms. This is because the ``atoms()`` function
returns a list of cursors, one per atom.

>>> cursor = mol.cursor()
>>> print(cursor.atoms())
Cursors( size=22
1: Cursor(atom, HH31:1)
2: Cursor(atom, CH3:2)
3: Cursor(atom, HH32:3)
4: Cursor(atom, HH33:4)
5: Cursor(atom, C:5)
...
18: Cursor(atom, H:18)
19: Cursor(atom, CH3:19)
20: Cursor(atom, HH31:20)
21: Cursor(atom, HH32:21)
22: Cursor(atom, HH33:22)
)
>>> for atom in cursor.atoms():
...     atom["color"] = atom["element"].color_name()
>>> mol = cursor.commit()
>>> print(mol.property("color"))
SireMol::AtomStringProperty( size=22
0: white
1: charcoal
2: white
3: white
4: charcoal
...
17: white
18: charcoal
19: white
20: white
21: white
)

.. note::

    Note how we have used the :func:`sire.mol.Element.color_name`
    function of :func:`sire.mol.Element` to get the color typically
    used to represent each atom in a molecular viewer.

This process of creating a cursor, then applying a change to every single
atom in the cursor, then commiting the changes back the molecule, is very common.
It is so common, that sire provides the ``apply`` function to enable
you to write this as a single line of code;

>>> mol = mol.cursor().atoms().apply(
...     lambda atom: atom.set("color", atom["element"].color_name())).commit()
>>> print(mol.property("color"))
SireMol::AtomStringProperty( size=22
0: white
1: charcoal
2: white
3: white
4: charcoal
...
17: white
18: charcoal
19: white
20: white
21: white
)

.. note::

    Note how we have to use the ``atom.set("color", ...)`` rather than
    ``atom["color"] = ...`` in the lambda expression. This is because
    assignments (using ``=``) are not supported in a Python lambda expression.

Searching by property
---------------------

You have :doc:`already seen <../part02/10_searching>` how to search for the more
standard properties, such as ``element``, ``mass`` and ``charge``.

You can also search for custom properties, such as the ``color`` property
we added above, using ``atom property``.

>>> print(mol["atom property color == charcoal"])
Selector<SireMol::Atom>( size=6
0:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
1:  Atom( C:5     [  18.48,    4.55,   14.35] )
2:  Atom( CA:9    [  16.54,    5.03,   15.81] )
3:  Atom( CB:11   [  16.05,    6.39,   15.26] )
4:  Atom( C:15    [  15.37,    4.19,   16.43] )
5:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
)

This supports any properties that are strings, numbers or boolean types.
All of the standard comparison operators (e.g. ``==``, ``>=``, ``!=`` etc.)
are supported.

For example, we could add a ``radius`` property based on each element's
covalent radius...

>>> cursor = mol.cursor()
>>> for atom in cursor.atoms():
...    atom["radius"] = atom["element"].covalent_radius().value()
>>> mol = cursor.commit()
>>> print(mol.property("radius"))
SireMol::AtomFloatProperty( size=22
0: 0.23
1: 0.68
2: 0.23
3: 0.23
4: 0.68
...
17: 0.23
18: 0.68
19: 0.23
20: 0.23
21: 0.23
)

.. note::

    Note how we have used the ``.value()`` function on the radius
    to get the raw value of the radius, without the units.
    We need to do this because we want to be able to search using
    the radius. Searching can only be performed with simple (numeric
    or boolean) properties.

or, using ``apply``, this could be written as

>>> mol = mol.cursor().atoms().apply(
...     lambda atom: atom.set("radius",
...                           atom["element"].covalent_radius().value()
...                          )).commit()

.. note::

    You can use either ``.cursor().atoms()`` or ``.atoms().cursor()`` - the
    order does not change the result.

We can now search for all atoms that have a radius that is less than ``0.5``.

>>> print(mol["atom property radius < 0.5"])
Selector<SireMol::Atom>( size=12
0:  Atom( HH31:1  [  18.45,    3.49,   12.44] )
1:  Atom( HH32:3  [  20.05,    3.63,   13.29] )
2:  Atom( HH33:4  [  18.80,    2.43,   13.73] )
3:  Atom( H:8     [  16.68,    3.62,   14.22] )
4:  Atom( HA:10   [  17.29,    5.15,   16.59] )
...
7:  Atom( HB3:14  [  15.24,    6.18,   14.55] )
8:  Atom( H:18    [  15.34,    5.45,   17.96] )
9:  Atom( HH31:20 [  14.35,    3.41,   19.15] )
10:  Atom( HH32:21 [  13.19,    4.59,   18.94] )
11:  Atom( HH33:22 [  13.21,    3.33,   17.69] )
)

Boolean properties are particularly useful, as these can be used to
mark atoms as matching particular criteria.

For example, we could set a property that is ``True`` for oxygen atoms using either

>>> cursor = mol.cursor()
>>> for atom in cursor.atoms("element O"):
...     atom["special"] = True
>>> mol = cursor.commit()

or

>>> mol = mol.cursor().atoms("element O").apply(
...             lambda atom: atom.set("special", True)).commit()

and can then use this property to search for those atoms.

>>> print(mol["atom property special == True"])
Selector<SireMol::Atom>( size=2
0:  Atom( O:6     [  19.19,    5.44,   14.76] )
1:  Atom( O:16    [  14.94,    3.17,   15.88] )
)

This search can be simplified to

>>> print(mol["atom property special"])
Selector<SireMol::Atom>( size=2
0:  Atom( O:6     [  19.19,    5.44,   14.76] )
1:  Atom( O:16    [  14.94,    3.17,   15.88] )
)

This is because an atom property search will return all of the atoms
that have a non-zero, non-empty or non-false value for the specified
property.

Deleting properties
-------------------

You can remove properties from the cursor in the same way that you
remove properties from a normal Python dictionary. You just ``del``
the key for the property you want to remove, or call the
``delete`` function of the Cursor, passing in the key.

For example, we can delete the ``radius`` property we created earlier
using

>>> cursor = mol.cursor()
>>> del cursor["radius"]
>>> mol = cursor.commit()
>>> print(mol.property("radius"))
KeyError: 'SireBase::missing_property: There is no property with
name "radius". Available properties are [ velocity, element, gb_radius_set,
bond, forcefield, gb_radii, color, angle, improper, gb_screening, intrascale,
LJ, coordinates, dihedral, treechain, connectivity, special, charge,
ambertype, parameters, atomtype, mass ].
(call sire.error.get_last_error_details() for more info)'

Or, alternatively, using the ``delete`` function,

>>> cursor = mol.cursor()
>>> cursor.delete("radius")
>>> mol = cursor.commit()

or, as one line,

>>> mol = mol.cursor().delete("radius").commit()

We can also remove the properties from individual atoms. Here, we will
remove the ``special`` property from the oxygen atoms

>>> cursor = mol.cursor()
>>> for atom in cursor.atoms("element O"):
...     del atom["special"]
>>> mol = cursor.commit()

or, alternatively, using the ``delete`` function,

>>> mol = mol.cursor().atoms("element O").delete("special").commit()

Deleting a property from an atom will reset it to a default-constructed
value. This is ``False`` (or ``0``) for boolean properties.

>>> print(mol["element O"].property("special"))
[0, 0]

While this is what you want for boolean properties, this may give
unexpected results for more complex properties. For example, deleting
the ``coordinates`` property from an atom will set its coordinates to
``(0, 0, 0)``...

>>> mol = mol.cursor().atoms("element O").delete("coordinates").commit()
>>> print(mol["element O"].property("coordinates"))
[( 0, 0, 0 ), ( 0, 0, 0 )]
