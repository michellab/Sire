"""
.. currentmodule:: sire.units

"""

from ..legacy import Units as _Units

from .. import use_new_api as _use_new_api
_use_new_api()

from ..legacy.Units import *


def _fix_generalunit():
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


if not hasattr(GeneralUnit, "approx_equal"):
    _fix_generalunit()


def switch_to_si():
    """
    Switch over to using SI units for output and default conversions

    This uses:

        mass:    gram (g)
        energy:  kilojoule (kJ)

    """
    gram.set_as_default("g")
    kilojoule.set_as_default("kJ")
    nanometer.set_as_default("nm")
    kJ_per_mol.set_as_default("kJ mol-1")
    mole.set_as_default("mol")
    picosecond.set_as_default("ps")


def switch_to_akma():
    gram.set_as_default("g")
    kcal.set_as_default("kcal")
    angstrom.set_as_default("Å")
    kcal_per_mol.set_as_default("kcal mol-1")
    mole.set_as_default("mol")
    picosecond.set_as_default("ps")
