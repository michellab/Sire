===========
Sire.System
===========

This module provides the classes that are used to collect and
manage molecular, forcefield and other data into complete
systems (held in :class:`Sire.System.System` objects).
The :class:`~Sire.System.System` is the complete molecular
system, with all information needed for simulation.

The moves in :mod:`Sire.Move` operate on :class:`~Sire.System.System`
objects. These are also what are returned by the molecule
parsers in :mod:`Sire.IO`.

As for much of Sire, you rarely need to instantiate or use these
classes directly. They are created and managed for you by
higher level parts of Sire.

.. toctree::
   :maxdepth: 3

   index_api_Sire_System
