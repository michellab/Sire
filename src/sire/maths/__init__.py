"""
.. currentmodule:: sire.maths


"""

from ..legacy import Maths as _Maths

from .. import use_new_api as _use_new_api
_use_new_api()

Triangle = _Maths.Triangle
Torsion = _Maths.Torsion

pi = _Maths.pi

from ._vector import Vector
from ._sphere import Sphere
