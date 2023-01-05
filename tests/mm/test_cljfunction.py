import pytest

import sire as sr


def test_cljfunction(ala_mols):
    mols = ala_mols

    # these are the two closest water molecules
    mol0 = mols[187]
    mol1 = mols[611]

    # This is how (at a low level) we can convert properties
    # into CLJAtoms (vector of coordinates, charges and LJ parameters)
    # We set the molecule ID to 1
    clj0 = sr.legacy.MM.CLJAtoms(
        mol0.property("coordinates").to_vector(),
        mol0.property("charge").to_vector(),
        mol0.property("LJ").to_vector(),
        1,
    )

    # Same for the other molecule, except with molecule ID 2
    clj1 = sr.legacy.MM.CLJAtoms(
        mol1.property("coordinates").to_vector(),
        mol1.property("charge").to_vector(),
        mol1.property("LJ").to_vector(),
        2,
    )

    # Now create a CLJFunction to calculate CLJ energies with
    # 15 A cutoff
    func = sr.legacy.MM.CLJShiftFunction(15 * sr.units.angstrom)

    nrgs = func.calculate(clj0, clj1)

    assert nrgs[0] == pytest.approx(-10.90987, 3)
    assert nrgs[1] == pytest.approx(7.1808714, 3)
