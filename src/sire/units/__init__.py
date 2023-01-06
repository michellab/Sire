__all__ = [
    "angle",
    "convert",
    "length",
    "set_energy_unit",
    "set_internal_units",
    "set_length_unit",
    "set_mass_unit",
    "set_quantity_unit",
    "set_si_units",
    "set_time_unit",
    "angstrom",
    "celsius",
    "degrees",
    "fahrenheit",
    "farad",
    "femtosecond",
    "g_per_mol",
    "kcal_per_mol",
    "kJ_per_mol",
    "mod_electron",
    "meter",
    "nanometer",
    "nanometer2",
    "nanosecond",
    "picometer",
    "picosecond",
    "radians",
]

from ..legacy import Units as _Units

from .. import use_new_api as _use_new_api

from ..legacy.Units import (
    angstrom,
    cal,
    celsius,
    convert,
    degrees,
    fahrenheit,
    farad,
    femtosecond,
    g_per_mol,
    gram,
    joule,
    mod_electron,
    meter,
    int_cal,
    int_kcal,
    kcal,
    kcal_per_mol,
    kilogram,
    kilojoule,
    kJ_per_mol,
    mole,
    nanometer,
    nanometer2,
    nanosecond,
    picometer,
    picosecond,
    radians,
    GeneralUnit,
)

_use_new_api()


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

    def __generalunit__getitem__(obj, key):
        return obj.get_component(key)

    def __generalunit__setitem__(obj, key, value):
        obj.set_component(key, value)
        return obj

    GeneralUnit.approx_equal = _generalunit_approx_equal
    GeneralUnit.to_default = _generalunit_to_default
    GeneralUnit.__setitem__ = __generalunit__setitem__
    GeneralUnit.__getitem__ = __generalunit__getitem__

    def __generalunit__bool__(obj):
        return not obj.is_zero()

    def __generalunit__float__(obj):
        if not obj.is_dimensionless():
            raise TypeError(
                f"You cannot convert the dimensioned value {obj} "
                "to a dimensionless floating point number."
            )

        return obj.value()

    def __generalunit__int__(obj):
        if not obj.is_dimensionless():
            raise TypeError(
                f"You cannot convert the dimensioned value {obj} "
                "to a dimensionless whole number integer."
            )

        return int(obj.value())

    GeneralUnit.__bool__ = __generalunit__bool__
    GeneralUnit.__float__ = __generalunit__float__
    GeneralUnit.__int__ = __generalunit__int__


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
            femtosecond: "fs",
        }

    try:
        return _names[unit]
    except Exception:
        pass

    return unit._to_cpp_type()


def _incompatible_units(a, b, typ: str):
    raise TypeError(
        f"Unit {a.unit_string()} is not a {typ}, and so is "
        f"incompatible with {b.unit_string()}."
    )


def set_quantity_unit(q, name: str = None):
    """Set the default quantity unit to be used for
    output and default conversions
    """

    if not q.has_same_units(mole):
        _incompatible_units(q, mole, "quantity")

    if name is None:
        name = _get_unit_name(q)

    q.set_as_default(name)

    e = kcal.get_default()
    ename = e.unit_string()

    (e / q).set_as_default(f"{ename} {name}-1")

    length = angstrom.get_default()
    lname = length.unit_string()

    (e / (q * length)).set_as_default(f"{ename} {name}-1 {lname}-1")
    (e / (q * length * length)).set_as_default(f"{ename} {name}-1 {lname}-2")

    g = gram.get_default()
    gname = g.unit_string()

    (g / q).set_as_default(f"{gname} {name}-1")


def set_energy_unit(energy, name: str = None):
    """Set the default energy unit to be used for
    output and default conversions
    """

    if not energy.has_same_units(kcal):
        _incompatible_units(energy, kcal, "energy")

    if name is None:
        name = _get_unit_name(energy)

    q = mole.get_default()
    qname = q.unit_string()

    energy.set_as_default(name)
    (energy / q).set_as_default(f"{name} {qname}-1")


def set_length_unit(length, name: str = None):
    """Set the default length unit to be used for
    output and default conversions
    """

    if not length.has_same_units(angstrom):
        _incompatible_units(length, angstrom, "length")

    if name is None:
        name = _get_unit_name(length)

    length.set_as_default(name)
    (length * length).set_as_default(f"{name}^2")
    (length * length * length).set_as_default(f"{name}^3")

    (1.0 / length).set_as_default(f"{name}-1")
    (1.0 / (length * length)).set_as_default(f"{name}-2")
    (1.0 / (length * length * length)).set_as_default(f"{name}-3")

    e = kcal_per_mol.get_default()
    ename = e.unit_string()

    (e / length).set_as_default(f"{ename} {name}-1")
    (e / (length * length)).set_as_default(f"{ename} {name}-2")


def set_mass_unit(mass, name: str = None):
    """Set the default mass unit to be used for
    output and default conversions
    """

    if not mass.has_same_units(gram):
        _incompatible_units(mass, gram, "mass")

    if name is None:
        name = _get_unit_name(mass)

    mass.set_as_default(name)

    q = mole.get_default()
    qname = q.unit_string()

    (mass / q).set_as_default(f"{name} {qname}-1")


def set_time_unit(time, name: str = None):
    """Set the default time unit to be used for
    output and default conversions
    """

    if not time.has_same_units(picosecond):
        _incompatible_units(time, picosecond, "time")

    if name is None:
        name = _get_unit_name(time)

    time.set_as_default(name)

    (1.0 / time).set_as_default(f"{name}-1")
    (1.0 / (time * time)).set_as_default(f"{name}-2")

    length = angstrom.get_default()
    lname = length.unit_string()

    (length / time).set_as_default(f"{lname} {name}-1")
    (length / (time * time)).set_as_default(f"{lname} {name}-2")


def set_si_units():
    """
    Switch over to using SI units for output and default conversions

    This uses:

        mass:     gram (g)
        energy:   kilojoule (kJ)
        length:   nanometer
        time:     picosecond
        quantity: mole
        charge:   mod_electron charges
    """
    set_quantity_unit(mole)
    set_energy_unit(kilojoule)
    set_mass_unit(gram)
    set_time_unit(picosecond)
    set_length_unit(nanometer)


def set_internal_units():
    """
    Switch over to using default (AKMA-style) units for output
    and default conversions.

    This uses:

        mass:     gram (g)
        energy:   kilocalorie (kcal)
        length:   angstrom
        time:     picosecond
        quantity: mole
        charge:   mod_electron charges
    """
    set_quantity_unit(mole)
    set_energy_unit(kcal)
    set_mass_unit(gram)
    set_time_unit(picosecond)
    set_length_unit(angstrom)


def length(length):
    """Convert the passed argument into a length. This will
    return `length` if it is already a length, or it will
    convert it into a length with default units if this
    is a float or something that can be converted to
    a float
    """
    try:
        return float(length) * angstrom.get_default()
    except Exception:
        pass

    if type(length) != type(angstrom):
        raise TypeError(
            f"The value '{length}' of type {type(length)} is "
            "not a type that is compatible with a GeneralUnit Length"
        )

    return length


def angle(angle):
    """Convert the passed argument into an angle. This will
    return 'angle' if it is already an angle, or will
    convert it into an angle with default units if this
    is a float or something that can be converted to
    a float
    """
    try:
        return float(angle) * degrees.get_default()
    except Exception:
        pass

    if type(angle) != type(degrees):
        raise TypeError(
            f"The value '{angle}' of type {type(angle)} is "
            "not a type that is compatible with a GeneralUnit angle"
        )

    return angle
