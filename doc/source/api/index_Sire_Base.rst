=========
Sire.Base
=========

This module implements base functionality that is used across Sire.
A key part of this is :class:`sire.base.Property`. Nearly all
objects in Sire derive from :class:`~sire.base.Property`. These can
be held as generic collections, e.g. a :class:`sire.mol.Molecule` is
really just a collection of properties, keyed as `coordinates`,
`charge` and `element`.

:func:`~sire.base.Property`
    The base class of most Sire object.

:func:`~sire.base.Properties`
    Used to hold collections of :class:`~sire.base.Property` objects,
    keyed by :class:`~sire.base.PropertyName`.

.. toctree::
   :maxdepth: 3

   index_api_Sire_Base
