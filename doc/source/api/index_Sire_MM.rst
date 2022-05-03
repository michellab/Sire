=======
Sire.MM
=======

This module implements the molecular mechanics forcefields. These
are forcefields that calculate molecular energies and forces
using classical molecular mechanics potentials.

There are a lot of classes in this module. Many of them are only
used internally as part of Sire's use as a rapid prototyping
engine.

The useful classes are those that implement complete molecular
mechanics forcefields. Sire does this by using forcefield classes
that represent different forcefield components:

:class:`~sire.mm.InterCLJFF`
    Calculates the intermolecular coulomb and Lennard Jones (CLJ)
    energy of all contained molecules.

:class:`~sire.mm.InterGroupCLJFF`
    Calculates the intermolecular coulomb and Lennard Jones (CLJ)
    energy between the two groups of molecules contained within.

:class:`~sire.mm.IntraCLJFF`
    Calculates the intramolecular coulomb and Lennard Jones (CLJ)
    energy within all contained molecules.

:class:`~sire.mm.IntraGroupCLJFF`
    Calculates the intramolecular coulomb and Lennard Jones (CLJ)
    energy between groups within a molecule.

:class:`~sire.mm.InternalFF`
    Calculates the intramolecular bond, angle and dihedral energy
    of contained molecules. These energies are calculated via
    algebraic expressions.

Each forcefield exposes its energy components as :mod:`sire.cas`
algebraic symbols. You can then assemble a total energy by combining
these symbols in whatever way you need.

Sire was originally written as a means to prototype new Monte Carlo
algorithms. As such, calculations of energies are highly optimised
for Monte Carlo. Calculations of changes in energy resulting from
small changes in molecular positions or properties are handled
in as efficient a manner as possible. Because of this, the forcefield
classes are less efficient at calculating forces or performing
molecular dynamics. We highly recommend `OpenMM <https://openmm.org>`__
as a great library for prototyping molecular dynamics simulations.
We do have, and continue to work on a bridge between Sire and OpenMM.

.. toctree::
   :maxdepth: 3

   index_api_Sire_MM
