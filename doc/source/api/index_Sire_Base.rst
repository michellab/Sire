=========
Sire.Base
=========

This module implements base functionality that is used across Sire.
A key part of this is :class:`Sire.Base.Property`. Nearly all
objects in Sire derive from :class:`~Sire.Base.Property`. These can
be held as generic collections, e.g. a :class:`Sire.Mol.Molecule` is
really just a collection of properties, keyed as `coordinates`,
`charge` and `element`.

:func:`~Sire.Base.Property`
    The base class of most Sire object.

:func:`~Sire.Base.Properties`
    Used to hold collections of :class:`~Sire.Base.Property` objects,
    keyed by :class:`~Sire.Base.PropertyName`.

.. toctree::
   :maxdepth: 1

   index_api_Sire_Base
