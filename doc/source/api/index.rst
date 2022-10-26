=============
Documentation
=============

Sire is composed of several modules:

* The :doc:`Sire <index_Sire>` top-level Python module. This has top-level
  functions that should make it easy to use Sire without knowing about
  or directly loading the other modules.
  By convention, we normally import `Sire` under
  the alias ``sr``.

There are also a number of submodules. These provide additional functionality,
going all the way up to letting you construct complete molecular simulation
programs.

* The :doc:`sire.mol <index_Sire_Mol>` module. This contains all of the
  classes and functions used to represent and manipulate molecules.

* The :doc:`sire.system <index_Sire_System>` module. This contains
  classes which represent systems (collections) of molecules. Systems
  can combine molecules with other data, e.g. forcefields for calculating
  energies, or additional properties such as the periodic box holding
  the molecules.

* The :doc:`sire.io <index_Sire_IO>` module. This provides lots of
  parsers to read and write molecular information to a range of
  different file formats (e.g. PDB, mol2, grotop etc.).

* The :doc:`sire.units <index_Sire_Units>` module. This contains
  classes and functions for representing the units of physical
  quantities.

* The :doc:`sire.cas <index_Sire_CAS>` module. This implements a
  complete Computer Algebra System, which is used for, e.g. giving
  the algebraic form of bond or dihedral functions.

* The :doc:`sire.move <index_Sire_Move>` module. This provides a
  collection of moves, which can be used to move atoms or molecules
  in a :class:`~sire.system.System` as part of a simulation.

* The :doc:`sire.analysis <index_Sire_Analysis>` module. This
  provides a set of analysis functionality that can be used
  to analyse running simulations (e.g. calculating FEP and TI
  free energies).

* The :doc:`sire.ff <index_Sire_FF>` module. This provides the base
  functionality for all forcefields (classes used to calculate the
  energy or forces on molecules).

* The :doc:`sire.mm <index_Sire_MM>` module. This provides code to calculate
  energies and forces according to several molecular mechanics forcefields.

* The :doc:`sire.squire <index_Sire_Squire>` module. This provides
  interfaces to various quantum chemical packages, so that you can
  calculate QM or QM/MM energies.

* Other useful modules, such as; :doc:`sire.maths <index_Sire_Maths>`
  (maths routines); :doc:`sire.base <index_Sire_Base>` (base utilities,
  including the :class:`~sire.base.Property` class);
  :doc:`sire.stream <index_Sire_Stream>` (functions relating to streaming/pickling
  of Sire objects); :doc:`sire.id <index_Sire_ID>` (functions used to
  implement Sire's identification system) and
  :doc:`sire.vol <index_Sire_Vol>` (classes that represent different
  simulation spaces, e.g. :class:`sire.vol.PeriodicBox`).

Detail
------

.. toctree::
   :maxdepth: 1

   index_Sire
   index_Sire_Analysis
   index_Sire_Base
   index_Sire_CAS
   index_Sire_Error
   index_Sire_FF
   index_Sire_ID
   index_Sire_IO
   index_Sire_Maths
   index_Sire_MM
   index_Sire_Mol
   index_Sire_Move
   index_Sire_Squire
   index_Sire_Stream
   index_Sire_System
   index_Sire_Units
   index_Sire_Vol
