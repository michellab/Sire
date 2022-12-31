
from .. import FF as _FF
from .. import CAS as _CAS
from .. import Mol as _Mol

# Import all of the classes and functions from the C++ library
from ._MM import *

# Now define some pure Python functions and classes that are part of
# this library...

# Next define all of the MM forcefield types so that the code can
# get them
def _createMMTypes():
    amberff = MMDetail(name = "amber::ff",
                       combining_rules = "arithmetic",
                       scale14elec = 1.0/1.2, scale14vdw = 0.5,
                       elecstyle = "coulomb", vdwstyle = "lj",
                       bondstyle = "harmonic", anglestyle = "harmonic",
                       dihedralstyle = "cosine")

    amberff99 = MMDetail(name = "amber::ff99",
                         combining_rules = "arithmetic",
                         scale14elec = 1.0/1.2, scale14vdw = 0.5,
                         elecstyle = "coulomb", vdwstyle = "lj",
                         bondstyle = "harmonic", anglestyle = "harmonic",
                         dihedralstyle = "cosine")

    amberff03 = MMDetail(name = "amber::ff03",
                         combining_rules = "arithmetic",
                         scale14elec = 1.0/1.2, scale14vdw = 0.5,
                         elecstyle = "coulomb", vdwstyle = "lj",
                         bondstyle = "harmonic", anglestyle = "harmonic",
                         dihedralstyle = "cosine")

    amberff12 = MMDetail(name = "amber::ff12",
                         combining_rules = "arithmetic",
                         scale14elec = 1.0/1.2, scale14vdw = 0.5,
                         elecstyle = "coulomb", vdwstyle = "lj",
                         bondstyle = "harmonic", anglestyle = "harmonic",
                         dihedralstyle = "cosine")

    amberff14 = MMDetail(name = "amber::ff14",
                         combining_rules = "arithmetic",
                         scale14elec = 1.0/1.2, scale14vdw = 0.5,
                         elecstyle = "coulomb", vdwstyle = "lj",
                         bondstyle = "harmonic", anglestyle = "harmonic",
                         dihedralstyle = "cosine")

    ambergaff = MMDetail(name = "amber::gaff",
                         combining_rules = "arithmetic",
                         scale14elec = 1.0/1.2, scale14vdw = 0.5,
                         elecstyle = "coulomb", vdwstyle = "lj",
                         bondstyle = "harmonic", anglestyle = "harmonic",
                         dihedralstyle = "cosine")

# initialise all of the MM forcefield types
_createMMTypes()

def getForceFields():
    """Return the names of the different MM forcefields that are recognised by
       this program"""
    return _FF_FFDetail.forcefields()

def getForceField(name):
    """Return the MM forcefield called 'name'"""
    return _FF.FFDetail.get(name)

###### PROPERTY KLUDGE FIX

__props = [ AtomLJs ]

for __prop in __props:
    _Mol._pvt_property_cludge_fix(__prop)
