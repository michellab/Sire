"""
.. currentmodule:: sire.legacy.Vol

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

from .. import Base as _Base
from .. import Maths as _Maths
from .. import Units as _Units

# Import all of the classes and functions from the C++ library
from ._Vol import *

# Now define some pure Python functions and classes that are part of
# this library...


