"""
.. currentmodule:: sire.legacy.FF

Classes
=======

.. autosummary::
    :toctree: generated/

    AtomPoint
    Center
    CenterOfGeometry
    CenterOfMass
    EnergyTable
    FF
    FF3D
    FFComponent
    FFID
    FFIdx
    FFMolGroup
    FFName
    FieldTable
    ForceTable
    G1FF
    G2FF
    GridFieldTable
    GridPotentialTable
    MolEnergyTable
    MolFieldTable
    MolForceTable
    MolPotentialTable
    Point
    PointRef
    Probe
    SingleComponent
    VectorPoint

"""

from .. import Mol as _Mol
from .. import CAS as _CAS

# Import all of the classes and functions from the C++ library
from ._FF import *

# Now define some pure Python functions and classes that are part of
# the module

