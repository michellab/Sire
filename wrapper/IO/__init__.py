
from .. import Mol as _Mol
from .. import System as _System
from .. import Units as _Units
from .. import MM as _MM
from .. import Qt as _Qt

# Import all of the classes and functions from the C++ library
from ._IO import *

# This is now a hack that is used to fix the AtomProperty
# wrapping issues on Mac. For some reason, some AtomProperty
# classes are exposed to that their type is wrong, meaning
# that you can't call any of their member functions.
# Editing a property of this type magically fixes everything...
def _fix_atomproperty_types():
    """This will attempt to fix the issues relating to wrapping of
       AtomProperty types. Essentially, after sire.legacy.IO has been loaded,
       we need to create and then extract an AtomProperty of each type.
       This must be done before any AtomProperties are created in Python.

       If I have missed an AtomProperty type then just add it below
       and assert that the resulting .nAtoms() function is called
       correctly.
    """
    from ..System import create_test_molecule
    from ..Mol import AtomIdx, ResIdx, ChainIdx, SegIdx
    from ..Units import g_per_mol, angstrom

    m = create_test_molecule()

    try:
        m = m.edit() \
                .atom(AtomIdx(0)) \
                .setProperty("a_m", 5 * g_per_mol) \
                .setProperty("a_s", "CA") \
                .setProperty("a_f", 5.0) \
                .setProperty("a_l", 5 * angstrom).molecule() \
                .residue(ResIdx(0)) \
                .setProperty("r_s", "A") \
                .setProperty("r_f", 5.0).molecule() \
                .chain(ChainIdx(0)) \
                .setProperty("c_s", "A") \
                .setProperty("c_f", 5.0).molecule() \
                .segment(SegIdx(0)) \
                .setProperty("s_s", "A") \
                .setProperty("s_f", 5.0).molecule().commit()
    except Exception:
        # try the new API
        m = m.edit() \
                .atom(AtomIdx(0)) \
                .set_property("a_m", 5 * g_per_mol) \
                .set_property("a_s", "CA") \
                .set_property("a_f", 5.0) \
                .set_property("a_l", 5 * angstrom).molecule() \
                .residue(ResIdx(0)) \
                .set_property("r_s", "A") \
                .set_property("r_f", 5.0).molecule() \
                .chain(ChainIdx(0)) \
                .set_property("c_s", "A") \
                .set_property("c_f", 5.0).molecule() \
                .segment(SegIdx(0)) \
                .set_property("s_s", "A") \
                .set_property("s_f", 5.0).molecule().commit()

    try:
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
    except Exception:
        # try the new API
        assert m.property("a_m").num_atoms() == 1
        assert m.property("a_s").nun_atoms() == 1
        assert m.property("a_f").num_atoms() == 1
        assert m.property("a_l").num_atoms() == 1
        assert m.property("r_s").num_residues() == 1
        assert m.property("r_f").num_residues() == 1
        assert m.property("c_s").num_chains() == 1
        assert m.property("c_f").num_chains() == 1
        assert m.property("s_s").num_segments() == 1
        assert m.property("s_f").num_segments() == 1

try:
    _fix_atomproperty_types()
except Exception as e:
    print(e)
    print("WARNING: AtomProperty classes aren't wrapped correctly.")
    print("WARNING: This will cause some problems when manipulating the ")
    print("WARNING: properties of atoms in Python.")

# Now define some pure Python functions and classes that are part of
# this library...
