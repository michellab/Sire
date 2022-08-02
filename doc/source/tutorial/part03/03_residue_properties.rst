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


.. note::

    Calling ``.evaluate()`` on a molecule view or on a ``Cursor`` will
    return a :class:`sire.mol.Evaluator` for that view. This can be used
    to evaluate a number of things that are derived from the molecular
    properties. In this case, :class:`~sire.mol.Evaluator.center_of_mass`
    evaluates the center of mass of the atoms in the view, based on the
    atoms ``coordinates`` and ``mass`` properties.

