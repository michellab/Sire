"""
.. currentmodule:: sire.units

"""

from ..legacy import Units as _Units

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Units import *


def _generalunit_approx_equal(u, v):
    if not hasattr(u, "what"):
        u = GeneralUnit(u)

    if not hasattr(v, "what"):
        v = GeneralUnit(v)

    if u.what().endswith("Property"):
        return _generalunit_approx_equal(u.value(), v)
    elif v.what().endswith("Property"):
        return _generalunit_approx_equal(u, v.value())

    # make sure that the units are the same
    if u.has_same_units(v):
        from ..search import approx_equal
        return approx_equal(u.value(), v.value())
    else:
        return False


GeneralUnit.approx_equal = _generalunit_approx_equal
