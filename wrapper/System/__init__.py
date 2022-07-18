"""
.. currentmodule:: sire.legacy.System

Classes
=======

.. autosummary::
    :toctree: generated/

    AngleComponent
    AssignerGroup
    ChargeConstraint
    CheckPoint
    CloseMols
    ComponentConstraint
    Constraint
    Constraints
    DihedralComponent
    DistanceComponent
    DoubleDistanceComponent
    EnergyMonitor
    FreeEnergyMonitor
    GeometryComponent
    IDAssigner
    IdentityConstraint
    MoleculeConstraint
    MonitorComponent
    MonitorComponents
    MonitorID
    MonitorIdx
    MonitorMonitor
    MonitorName
    MonitorProperty
    PerturbationConstraint
    PolariseCharges
    PolariseChargesFF
    PropertyConstraint
    SpaceWrapper
    SysID
    SysIdx
    SysName
    System
    SystemMonitor
    SystemMonitors
    TripleDistanceComponent
    VolMapMonitor
    WindowedComponent

Functions
=========

.. autosummary::
    :toctree: generated/

    create_test_molecule

"""

from .. import FF as _FF
from .. import Base as _Base

from ._System import *

System.__setProperty__ = System.setProperty
System.setProperty = _Base.__set_property__
