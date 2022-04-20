"""
.. currentmodule:: Sire.MM

Classes
=======

.. autosummary::
    :toctree: generated/

    AmberAngle
    AmberBond
    AmberDihedral
    AmberDihPart
    AmberNB14
    AngleComponent
    AngleParameterName
    AngleRestraint
    AngleSymbols
    AtomFunction
    AtomFunctions
    AtomLJs
    BendBendComponent
    BendBendParameterName
    BendBendSymbols
    BondComponent
    BondParameterName
    BondSymbols
    ChargeParameterName
    ChargeParameterName3D
    CHARMMSwitchingFunction
    CLJ14Group
    CLJAtom
    CLJAtoms
    CLJBox
    CLJBoxDistance
    CLJBoxes
    CLJBoxIndex
    CLJCalculator
    CLJComponent
    CLJCutoffFunction
    CLJDelta
    CLJExtractor
    CLJFunction
    CLJGrid
    CLJIntraFunction
    CLJIntraRFFunction
    CLJIntraShiftFunction
    CLJNBPairs
    CLJParameterNames
    CLJParameterNames3D
    CLJProbe
    CLJRFFunction
    CLJScaleFactor
    CLJShiftFunction
    CLJSoftFunction
    CLJSoftIntraFunction
    CLJSoftIntraRFFunction
    CLJSoftIntraShiftFunction
    CLJSoftRFFunction
    CLJSoftShiftFunction
    CLJWorkspace
    CoulombComponent
    CoulombNBPairs
    CoulombProbe
    CoulombScaleFactor
    DihedralComponent
    DihedralParameterName
    DihedralRestraint
    DihedralSymbols
    DistanceRestraint
    DoubleDistanceRestraint
    FourAtomFunction
    FourAtomFunctions
    FourAtomPerturbation
    GridFF
    GridFF2
    GromacsAngle
    GromacsAtomType
    GromacsBond
    GromacsDihedral
    GroupInternalParameters
    HarmonicSwitchingFunction
    ImproperComponent
    ImproperParameterName
    ImproperSymbols
    InterCLJFF
    InterCoulombFF
    InterFF
    InterGroupCLJFF
    InterGroupCoulombFF
    InterGroupFF
    InterGroupLJFF
    InterGroupSoftCLJFF
    InterLJFF
    InternalComponent
    InternalFF
    InternalParameterNames
    InternalParameters
    InternalParameters3D
    InternalPerturbation
    InternalSymbols
    InterSoftCLJFF
    Intra14Component
    Intra14CoulombComponent
    Intra14LJComponent
    IntraCLJFF
    IntraCoulombFF
    IntraFF
    IntraGroupCLJFF
    IntraGroupCoulombFF
    IntraGroupFF
    IntraGroupLJFF
    IntraGroupSoftCLJFF
    IntraLJFF
    IntraSoftCLJFF
    LJComponent
    LJNBPairs
    LJParameter
    LJParameterName
    LJParameterName3D
    LJPerturbation
    LJProbe
    LJScaleFactor
    MultiCLJComponent
    NoCutoff
    Restraint
    Restraint3D
    RestraintComponent
    RestraintFF
    ScaledChargeParameterNames3D
    ScaledCLJParameterNames3D
    ScaledLJParameterNames3D
    SoftCLJComponent
    StretchBendComponent
    StretchBendParameterName
    StretchBendSymbols
    StretchBendTorsionComponent
    StretchBendTorsionParameterName
    StretchBendTorsionSymbols
    StretchStretchComponent
    StretchStretchParameterName
    StretchStretchSymbols
    SwitchingFunction
    TestFF
    ThreeAtomFunction
    ThreeAtomFunctions
    ThreeAtomPerturbation
    TripleDistanceRestraint
    TwoAtomFunction
    TwoAtomFunctions

"""

import Sire.FF
import Sire.CAS

# Import all of the classes and functions from the C++ library
from Sire.MM._MM import *

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
    return Sire.FF.FFDetail.forcefields()

def getForceField(name):
    """Return the MM forcefield called 'name'"""
    return Sire.FF.FFDetail.get(name)

###### PROPERTY KLUDGE FIX

__props = [ AtomLJs ]

for __prop in __props:
    Sire.Mol._pvt_property_cludge_fix(__prop)
