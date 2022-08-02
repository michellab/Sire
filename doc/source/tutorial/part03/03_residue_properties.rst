=====================================
Residue, Chain and Segment Properties
=====================================

Residues, chains and segments can also have their own properties.

For example, we will now create a new residue property that holds the
radius of each residue. First, let's load the molecules again.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.crd"]))
>>> mol = mols[0]

Now, we will use a cursor to add a ``residue_center`` property to each
residue.

>>> cursor = mol.cursor()
>>> for residue in cursor.residues():
...     center = residue.evaluate().center_of_mass()
...     residue["residue_center"] = center
>>> mol = cursor.commit()
>>> print(mol.property("residue_center"))
SireMol::ResPropertyProperty( size=3
0: ( 18.9264, 4.47803, 14.1498 )
1: ( 16.0195, 4.60992, 15.5812 )
2: ( 14.3873, 4.27636, 18.0041 )
)

.. note::

    Calling ``.evaluate()`` on a molecule view or on a ``Cursor`` will
    return a :class:`sire.mol.Evaluator` for that view. This can be used
    to evaluate a number of things that are derived from the molecular
    properties. In this case, :class:`~sire.mol.Evaluator.center_of_mass`
    evaluates the center of mass of the atoms in the view, based on the
    atoms ``coordinates`` and ``mass`` properties.

We will also add a radius for each residue. We will do this by calling
:func:`sire.mol.Evaluator.bounding_sphere` to get a
:class:`sire.maths.Sphere` object that represents the bounding sphere
for each residue. We will then get the radius of this sphere by calling
:func:`sire.maths.Sphere.radius`.

We could do this either via a loop...

>>> cursor = mol.cursor()
>>> for residue in cursor.residues():
...    residue["residue_radius"] = residue.evaluate().bounding_sphere().radius()
>>> mol = cursor.commit()
>>> print(mol.property("residue_radius"))
SireMol::ResFloatProperty( size=3
0: 2.06211
1: 2.53194
2: 1.69918
)

...or via ``apply`` and a lambda expression...

>>> mol = mol.residues().cursor().apply(
...         lambda res: res.set("residue_radius",
...                             res.evaluate().bounding_sphere().radius())
...  ).commit()
>>> print(mol.property("residue_radius"))
SireMol::ResFloatProperty( size=3
0: 2.06211
1: 2.53194
2: 1.69918
)


