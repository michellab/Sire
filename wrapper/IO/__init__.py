"""
.. currentmodule:: Sire.IO

Classes
=======

.. autosummary::
    :toctree: generated/

    Amber
    AmberPrm
    AmberRst
    AmberRst7
    CharmmPSF
    Cube
    Gro87
    Mol2
    MoleculeParser
    PDB
    PDB2
    SDF
    Tinker
    TrajectoryMonitor
    ZmatrixMaker

Functions
=========

.. autosummary::
    :toctree: generated/

    isAmberWater
    isGromacsWater
    isWater
    renumberConstituents
    repartitionHydrogenMass
    setAmberWater
    setGromacsWater
    updateAndPreserveOrder
    updateCoordinatesAndVelocities

"""
import Sire.Mol
import Sire.System
import Sire.Units
import Sire.MM
import Sire.Qt

# Import all of the classes and functions from the C++ library
from Sire.IO._IO import *

# This is now a hack that is used to fix the AtomProperty
# wrapping issues on Mac. For some reason, some AtomProperty
# classes are exposed to that their type is wrong, meaning
# that you can't call any of their member functions.
# Editing a property of this type magically fixes everything...
def _fix_atomproperty_types():
    """This will attempt to fix the issues relating to wrapping of
       AtomProperty types. Essentially, after Sire.IO has been loaded,
       we need to create and then extract an AtomProperty of each type.
       This must be done before any AtomProperties are created in Python.

       If I have missed an AtomProperty type then just add it below
       and assert that the resulting .nAtoms() function is called
       correctly.
    """
    m = Sire.System.create_test_molecule()

    m = m.edit() \
            .atom(Sire.Mol.AtomIdx(0)) \
            .setProperty("a_m", 5 * Sire.Units.g_per_mol) \
            .setProperty("a_s", "CA") \
            .setProperty("a_f", 5.0) \
            .setProperty("a_l", 5 * Sire.Units.angstrom).molecule() \
            .residue(Sire.Mol.ResIdx(0)) \
            .setProperty("r_s", "A") \
            .setProperty("r_f", 5.0).molecule() \
            .chain(Sire.Mol.ChainIdx(0)) \
            .setProperty("c_s", "A") \
            .setProperty("c_f", 5.0).molecule() \
            .segment(Sire.Mol.SegIdx(0)) \
            .setProperty("s_s", "A") \
            .setProperty("s_f", 5.0).molecule().commit()

    assert m.property("a_m").nAtoms() == 1
    assert m.property("a_s").nAtoms() == 1
    assert m.property("a_f").nAtoms() == 1
    assert m.property("a_l").nAtoms() == 1
    assert m.property("r_s").nResidues() == 1
    assert m.property("r_f").nResidues() == 1
    assert m.property("c_s").nChains() == 1
    assert m.property("c_f").nChains() == 1
    assert m.property("s_s").nSegments() == 1
    assert m.property("s_f").nSegments() == 1

try:
    _fix_atomproperty_types()
except Exception as e:
    print(e)
    print("WARNING: AtomProperty classes aren't wrapped correctly.")
    print("WARNING: This will cause some problems when manipulating the ")
    print("WARNING: properties of atoms in Python.")

# Now define some pure Python functions and classes that are part of
# this library...

