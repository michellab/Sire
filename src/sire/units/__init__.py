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

    def _generalunit_to_default(obj):
        """Return a floating point value that represents
           this value in default units for this dimension

           Example
           -------

           >>> import sire as sr
           >>> l = 5 * sr.units.angstrom
           >>> sr.units.set_length_unit(sr.units.picometer)
           >>> print(l.to_default())
               500.0
        """
        return obj.value() / obj.get_default().value()


    GeneralUnit.approx_equal = _generalunit_approx_equal
    GeneralUnit.to_default = _generalunit_to_default


if not hasattr(GeneralUnit, "approx_equal"):
    _fix_generalunit()


_names = None

def _get_unit_name(unit):

    global _names

    if _names is None:
        _names = {
        mole: "mol",
        cal: "cal",
        joule: "joule",
        int_cal: "int_cal",
        kcal: "kcal",
        kilojoule: "kJ",
        int_kcal: "int_kcal",
        angstrom: "Å",
        picometer: "pm",
        nanometer: "nm",
        gram: "g",
        kilogram: "kg",
        picosecond: "ps",
        nanosecond: "ns",
        femtosecond: "fs"
        }

    try:
        return _names[unit]
    except Exception:
        pass

    return unit._to_cpp_type()


def _incompatible_units(a, b, typ:str):
    raise TypeError(
        f"Unit {a.unit_string()} is not a {typ}, and so is "
        f"incompatible with {b.unit_string()}."
    )


def set_quantity_unit(q, name: str=None):

    if not q.has_same_units(mole):
        _incompatible_units(q, mole, "quantity")

    if name is None:
        name = _get_unit_name(q)

    q.set_as_default(name)

    e = kcal.get_default()
    ename = e.unit_string()

    (e / q).set_as_default(f"{ename} {name}-1")

    l = angstrom.get_default()
    lname = l.unit_string()

    (e / (q*l)).set_as_default(f"{ename} {name}-1 {lname}-1")
    (e / (q*l*l)).set_as_default(f"{ename} {name}-1 {lname}-2")

    g = gram.get_default()
    gname = g.unit_string()

    (g / q).set_as_default(f"{gname} {name}-1")


def set_energy_unit(energy, name: str=None):

    if not energy.has_same_units(kcal):
        _incompatible_units(energy, kcal, "energy")

    if name is None:
        name = _get_unit_name(energy)

    q = mole.get_default()
    qname = q.unit_string()

    energy.set_as_default(name)
    (energy / q).set_as_default(f"{name} {qname}-1")


def set_length_unit(length, name: str=None):

    if not length.has_same_units(angstrom):
        _incompatible_units(length, angstrom, "length")

    if name is None:
        name = _get_unit_name(length)

    length.set_as_default(name)
    (length*length).set_as_default(f"{name}^2")
    (length*length*length).set_as_default(f"{name}^3")

    (1.0 / length).set_as_default(f"{name}-1")
    (1.0 / (length*length)).set_as_default(f"{name}-2")
    (1.0 / (length*length*length)).set_as_default(f"{name}-3")

    e = kcal_per_mol.get_default()
    ename = e.unit_string()

    (e / length).set_as_default(f"{ename} {name}-1")
    (e / (length*length)).set_as_default(f"{ename} {name}-2")


def set_mass_unit(mass, name:str = None):

    if not mass.has_same_units(gram):
        _incompatible_units(mass, gram, "mass")

    if name is None:
        name = _get_unit_name(mass)

    mass.set_as_default(name)

    q = mole.get_default()
    qname = q.unit_string()

    (mass / q).set_as_default(f"{name} {qname}-1")


def set_time_unit(time, name:str = None):

    if not time.has_same_units(picosecond):
        _incompatible_units(time, picosecond, "time")

    if name is None:
        name = _get_unit_name(time)

    time.set_as_default(name)

    (1.0 / time).set_as_default(f"{name}-1")
    (1.0 / (time*time)).set_as_default(f"{name}-2")

    l = angstrom.get_default()
    lname = l.unit_string()

    (l / time).set_as_default(f"{lname} {name}-1")
    (l / (time*time)).set_as_default(f"{lname} {name}-2")


def set_si_units():
    """
    Switch over to using SI units for output and default conversions

    This uses:

        mass:    gram (g)
        energy:  kilojoule (kJ)

    """
    set_quantity_unit(mole)
    set_energy_unit(kilojoule)
    set_mass_unit(gram)
    set_time_unit(picosecond)
    set_length_unit(nanometer)


def set_internal_units():
    set_quantity_unit(mole)
    set_energy_unit(kcal)
    set_mass_unit(gram)
    set_time_unit(picosecond)
    set_length_unit(angstrom)
