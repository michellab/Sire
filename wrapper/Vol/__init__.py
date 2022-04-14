"""
.. currentmodule:: Sire.Vol

Classes
=======

.. autosummary::
    :toctree: generated/

    AABox
    BoxPatching
    Cartesian
    CombinedSpace
    CombineSpaces
    CoordGroup
    CoordGroupArray
    CoordGroupArrayArray
    CoordGroupBase
    CoordGroupEditor
    Grid
    GridIndex
    GridInfo
    Patching
    PeriodicBox
    RegularGrid
    Space
    TriclinicBox

"""

import Sire.Base
import Sire.Maths
import Sire.Units

# Import all of the classes and functions from the C++ library
from Sire.Vol._Vol import *

__all__ = [ "AABox", "BoxPatching", "Cartesian", "CombinedSpace",
            "CombineSpaces", "CoordGroup", "CoordGroupArray", "CoordGroupArrayArray",
            "CoordGroupBase", "CoordGroupEditor", "Grid", "GridIndex",
            "GridInfo", "Patching", "PeriodicBox", "RegularGrid",
            "Space", "TriclinicBox" ]

# Now define some pure Python functions and classes that are part of
# this library...


