========
Sire.Mol
========

This module provides all of the molecular classes in Sire. The molecular
model in Sire is based on :class:`sire.base.Property` combined with
:class:`sire.mol.MoleculeInfo`.

A molecule is simply a collection of arbitrary properties, which are
indexed via the identities held in the attached :class:`~sire.mol.MoleculeInfo`.

This core data type (called ``MoleculeData`` in C++) is not exposed
directly. Instead, you interact with this via a series of types of
classes.

.. note::

    Note that you normally don't need to use the classes in this module
    directly. Instead, the classes will be created and used implicitly
    as you are calling higher-level functions.


View Classes
------------

You can view the data using one of the view classes. These include
:class:`~sire.mol.Atom`, :class:`~sire.mol.Residue`, :class:`~sire.mol.Chain`, :class:`~sire.mol.Segment` and :class:`~sire.mol.Molecule`, which provide views of atoms, residues,
chains, segments or whole molecules.

There are other view classes, e.g. :class:`~sire.mol.PartialMolecule`,
which can view an arbitrary selection of atoms. The atoms selected
are held in the :class:`~sire.mol.AtomSelection` class.

Edit Classes
------------

You can edit the data in the molecule by using one of the edit classes.
These are created by calling the ``.edit()`` function on a view.

This edits a copy of the molecule. Any edits have to be committed
back to the molecule, e.g.

.. code-block:: python

    >>> mol = mol.edit().setProperty("name", "My Molecule").commit()
    >>> print(mol.property("name"))
    "My Molecule"

We also provide a simpler interface to editing molecules via
the :class:`sire.Cursor` interface. This creates a cursor that automatically
handles the change to the underlying edit classes, e.g.

.. code-block:: python

    >>> c = mol.cursor()
    >>> c["name"] = "My Molecule"
    >>> print(c["name"])
    "My Molecule"

    >>> mol = c.commit()
    >>> print(mol.property("name"))
    "My Molecule"

Mover Classes
-------------

You can move a molecule (or part of a molecule) by using one of the move
classes. These work similarly to the edit classes, but instead have
functions to translate and rotate molecules (or parts of molecules).

.. code-block:: python

    >>> mol = mol.move().translate([1.0, 0.0, 2.0]).commit()

Evaluator Classes
-----------------

These classes perform calculations on the data held in a molecule
(or the part of the molecule being viewed). Examples include
calculating the total mass, charge or center of geometry.

.. code-block:: python

    >>> mass = mol.evaluate().mass()

Selector Classes
----------------

These provide views that hold selections within molecules. These selections
could be a set of Atoms, or Residues etc. You can index within molecules
using the ID classes, e.g. :class:`~sire.mol.AtomName` to locate by
atom name. Sire will do its best to work out what you mean if you just
pass in a string (likely a name) or a number (likely an index).

.. code-block:: python

    >>> mass = mol.atoms("CA").evaluate().mass()

You can also search for selections using the search function, e.g.

.. code-block:: python

    >>> mass = mol.search("atomname CA and resname ALA").evaluate().mass()

.. toctree::
   :maxdepth: 3

   index_api_Sire_Mol
