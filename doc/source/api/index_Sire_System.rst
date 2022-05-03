===========
Sire.System
===========

This module provides the classes that are used to collect and
manage molecular, forcefield and other data into complete
systems (held in :class:`sire.system.System` objects).
The :class:`~sire.system.System` is the complete molecular
system, with all information needed for simulation.

The moves in :mod:`sire.move` operate on :class:`~sire.system.System`
objects. These are also what are returned by the molecule
parsers in :mod:`sire.io`.

As for much of ``sire``, you rarely need to instantiate or use these
classes directly. They are created and managed for you by
higher level parts of ``sire``.

.. toctree::
   :maxdepth: 3

   index_api_Sire_System
