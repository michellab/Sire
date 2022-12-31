=============
Documentation
=============

:mod:`sire` is implemented as a primary module, with a collection
of submodules.

You import :mod:`sire` using;

>>> import sire as sr

This imports the main module, which is the primary gateway into all of
the functionality of :mod:`sire`.

This is designed so that you (mostly) only need to call functions
from this main module. These functions
will create objects and call functions from the various sub-modules.
You will rarely need to call functions or create classes from the
sub-modules directly yourself.

For example, :func:`sire.load` will use the file parsers defined
in :mod:`sire.io` to load molecules (each represented as a
:class:`~sire.mol.Molecule`
object from :mod:`sire.mol`), and will return the result as a
:class:`~sire.system.System` (defined in :mod:`sire.system`).

Key sub-modules are:

* :mod:`sire.mol` - defines all of the objects that are used to represent
  and manipulate molecules (and views of molecules).
* :mod:`sire.search` - provides the power behind the
  :doc:`search functionality <../cheatsheet/search>`.
* :mod:`sire.units` - implements a complete set of units so that all
  physical quantities are represented by proper units / dimensions.
* :mod:`sire.maths` - provides a collection of maths functions and classes
  (e.g. points in space, represented via :class:`sire.maths.Vector`, or
  spheres, represented via :class:`sire.maths.Sphere`).
* :mod:`sire.cas` - this is a complete computer algebra system (CAS), used
  to flexibly define interaction functions for calculating energies and
  forces between atoms and molecules.
* :mod:`sire.io` - all of the input/output classes used to read and write
  molecular data to and from files.
* :mod:`sire.mm` - a set of molecular mechanics (MM) forcefields that enables
  :mod:`sire` to rapidly calculate energies and forces between atoms and
  molecules.
* :mod:`sire.base` - a set of base classes and functions that includes a
  property system with associated containers for flexibly holding and
  manipulating objects of different types, and mapping them between
  the C++ and Python layers of :mod:`sire`.
* :mod:`sire.stream` - provides functions that support streaming (pickling)
  of all :mod:`sire` objects to an from a cross-platform, versioned binary
  format. This allows all C++ objects to be pickled via standard
  Python pickle.
* :mod:`sire.vol` - provides different spaces (volumes) that can be used
  to calculate distances between points with different boundary conditions.
* :mod:`sire.utils` - a miscellaneous collection of utilities that
  are Python-only, and don't fit neatly into any of the other sub-modules.

Top-level documentation
=======================

.. toctree::
   :maxdepth: 3

   sire

Sub-module documentation
========================

.. toctree::
   :maxdepth: 2

   sire_modules
