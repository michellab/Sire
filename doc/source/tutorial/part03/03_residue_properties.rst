=====================================
Residue, Chain and Segment Properties
=====================================

Residues, chains and segments can also have their own properties.

Residue properties
------------------

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
:func:`sire.mol.Evaluator.radius` and passing in the center calculated above.

We could do this either via a loop...

>>> cursor = mol.cursor()
>>> for residue in cursor.residues():
...    residue["residue_radius"] = residue.evaluate().radius(residue["residue_center"])
>>> mol = cursor.commit()
>>> print(mol.property("residue_radius"))

...or via ``apply`` and a lambda expression...

>>> mol = mol.residues().cursor().apply(
...         lambda res: res.set("residue_radius",
...                             res.evaluate().radius(res["residue_center"]))
...  ).commit()
>>> print(mol.property("residue_radius"))

Residue properties are a lot like atom properties. They behave like a list,
with one value per residue (in residue index order). Default-constructed
values are used for residues that don't have a property set.

Chain properties
----------------

Chain properties work in the same way as residue properties. For example,
load up protein ``7SA1``.

>>> mols = sr.load("7SA1")
>>> mol = mols[0]

We will loop over each chain, adding in a ``sphere`` property that contains
a bounding sphere for the entire chain.

We could do this via a loop...

>>> cursor = mol.cursor()
>>> for chain in cursor.chains():
...     chain["sphere"] = chain.evaluate().bounding_sphere()
>>> mol = cursor.commit()
>>> print(mol.property("sphere"))
SireMol::ChainPropertyProperty( size=4
0: Sphere( center() == ( -44.5915, 10.0365, 16.132 ), radius == 40.2531 )
1: Sphere( center() == ( -31.0725, 27.279, 3.387 ), radius == 65.6821 )
2: Sphere( center() == ( -5.827, 22.1885, 49.5585 ), radius == 40.2047 )
3: Sphere( center() == ( 6.3375, 5.929, 31.511 ), radius == 65.5556 )
)

...or via an ``apply``

>>> mol = mol.chains().cursor().apply(
...        lambda chain: chain.set("sphere",
...                                chain.evaluate().bounding_sphere())
...       ).commit()
>>> print(mol.property("sphere"))
SireMol::ChainPropertyProperty( size=4
0: Sphere( center() == ( -44.5915, 10.0365, 16.132 ), radius == 40.2531 )
1: Sphere( center() == ( -31.0725, 27.279, 3.387 ), radius == 65.6821 )
2: Sphere( center() == ( -5.827, 22.1885, 49.5585 ), radius == 40.2047 )
3: Sphere( center() == ( 6.3375, 5.929, 31.511 ), radius == 65.5556 )
)

Chain properties behave identically to residue properties, i.e. like
a list with one value per chain, in chain index order. Default-constructed
values are used for chains that don't have a property set.

Segment properties
------------------

Segment properties behave exactly like residue and chain properties. To
demonstrate this, load up ``alanin`` again...

>>> mols = sr.load(sr.expand(sr.tutorial_url, "alanin.psf"))
>>> mol = mols[0]

We will now add a property to the segments that is the number of atoms
in that segment. Again, we can do this via a loop...

>>> cursor = mol.cursor()
>>> for segment in cursor.segments():
...     segment["seg_atom_count"] = segment.view().num_atoms()
>>> mol = cursor.commit()
>>> print(mol.property("seg_atom_count"))
SireMol::SegIntProperty( size=1
0: 66
)

...or via an apply

>>> mol = mol.segments().cursor().apply(
...     lambda seg: seg.set("seg_atom_count", seg.view().num_atoms())
...   ).commit()
>>> print(mol.property("seg_atom_count"))
SireMol::SegIntProperty( size=1
0: 66
)

.. note::

    Note how we use ``cursor.view()`` to gain access to the view on which
    the cursor is operating. This is useful to call functions on that
    view, e.g. ``segment.num_atoms()`` in this case.

Segment properties behave identically to chain and residue properties, i.e.
like a list with one value per segment, in segment index order.
Default-constructed values are used for segments that don't have a
property set.
