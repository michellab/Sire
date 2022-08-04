
import pytest

from sire.units import angstrom, kcal_per_mol, GeneralUnit
from sire.base import Properties

def test_generalunit():
    l = 5 * angstrom

    g = GeneralUnit(l)

    assert g == l

    assert 2 * g == 2 * l

    assert g + l == l + g

    assert l / l == g / g

    assert g * g == l * l

    assert (g / g).is_dimensionless()

    assert g.to(angstrom) == l.value()

    with pytest.raises(UserWarning):
        g.to(kcal_per_mol)

    assert g.value() == l.value()
    assert g.units() == angstrom

    assert g.approx_equal(g)
    assert g.approx_equal(1.0000001 * g)
    assert g.approx_equal(0.9999999 * g)
    assert not g.approx_equal(1.1 * g)
    assert not g.approx_equal(0.9 * g)

    e = 4.4 * kcal_per_mol

    g = GeneralUnit(e)

    assert g == e
    assert 2 * g == 2 * e
    assert g + e == e + g
    assert e / e == g / g
    assert (g / g).is_dimensionless()

    assert g.to(kcal_per_mol) == e.value()

    with pytest.raises(UserWarning):
        g.to(angstrom)

    assert g.value() == e.value()
    assert g.units() == kcal_per_mol

    assert g.approx_equal(g)
    assert g.approx_equal(1.0000001 * g)
    assert g.approx_equal(0.9999999 * g)
    assert not g.approx_equal(1.1 * g)
    assert not g.approx_equal(0.9 * g)


def test_generalunitproperty():
    l = 5 * angstrom

    g = GeneralUnit(l)

    p = Properties()
    p["test"] = l

    h = p["test"]

    assert g == h

    assert g + g == h + h


def test_generalunit_zero():
    from sire.legacy.MM import CLJShiftFunction
    from sire.units import angstrom

    ff = CLJShiftFunction(100*angstrom, 150*angstrom)

    # this raises an exception is zero is not handled correctly
    ff = CLJShiftFunction(100*angstrom, 0*angstrom)


if __name__ == "__main__":
    test_generalunit()
    test_generalunitproperty()