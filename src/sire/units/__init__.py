"""
.. currentmodule:: sire.units

"""

from ..legacy import Units as _Units

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Units import *


from ..base import _add_property_operators
_add_property_operators(GeneralUnitProperty)
