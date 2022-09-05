
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


def test_generalunit_components():
    v = 0 * kcal_per_mol

    v.add_component("bond", 5 * kcal_per_mol)

    assert v == 5 * kcal_per_mol

    # note that this test validates that we have the
    # "principle of least surprise" equality that does
    # not compare the components of the energies
    assert v.get_component("bond") == 5 * kcal_per_mol

    with pytest.raises(UserWarning):
        v.set_component("angle", 3 * angstrom)

    v.add_component("angle", 3 * kcal_per_mol)

    assert v == 8 * kcal_per_mol

    assert v.get_component("bond") == 5 * kcal_per_mol
    assert v.get_component("angle") == 3 * kcal_per_mol

    v.add_component("bond", 10*kcal_per_mol)

    assert v == 18 * kcal_per_mol

    assert v.get_component("bond") == 15 * kcal_per_mol
    assert v.get_component("angle") == 3 * kcal_per_mol

    v.subtract_component("bond", 10*kcal_per_mol)

    assert v == 8 * kcal_per_mol

    assert v.get_component("bond") == 5 * kcal_per_mol
    assert v.get_component("angle") == 3 * kcal_per_mol

    v.set_component("bond", 2*kcal_per_mol)

    assert v == 5 * kcal_per_mol

    assert v.get_component("bond") == 2 * kcal_per_mol
    assert v.get_component("angle") == 3 * kcal_per_mol

    v.add_component("bond", -2*kcal_per_mol)

    assert v == 3 * kcal_per_mol

    assert v.get_component("bond") == 0
    assert v.get_component("angle") == 3 * kcal_per_mol

    assert v.get_component("dihedral") == 0

    v.set_component("angle", 0)

    assert v.get_component("bond") == 0
    assert v.get_component("angle") == 0
    assert v.get_component("dihedral") == 0

    assert v == 0

    v.set_component("bond", 5 * angstrom)

    assert v == 5 * angstrom

    assert v.get_component("bond") == 5 * angstrom

    with pytest.raises(UserWarning):
        v.set_component("angle", 3 * kcal_per_mol)

    v = GeneralUnit()

    v["bond"] = 5 * angstrom
    v["angle"] = 3 * angstrom

    assert v == 8 * angstrom

    assert v["bond"] == 5 * angstrom
    assert v["angle"] == 3 * angstrom

    v["bond"] += 2 * angstrom

    assert v["bond"] == 7 * angstrom
    assert v == 10 * angstrom

    v += v

    assert v == 20 * angstrom
    assert v["bond"] == 14 * angstrom
    assert v["angle"] == 6 * angstrom

    w = v * 5

    assert w == 100 * angstrom
    assert w["bond"] == 70 * angstrom
    assert w["angle"] == 30 * angstrom

    v *= 5

    assert v == 100 * angstrom
    assert v["bond"] == 70 * angstrom
    assert v["angle"] == 30 * angstrom

    assert v == w
    assert v.components() == w.components()
