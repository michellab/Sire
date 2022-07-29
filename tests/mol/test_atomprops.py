

import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_atomcoords(ala_mols):
    mols = ala_mols

    mol = mols[0]

    coords = mol.property("coordinates")

    assert len(coords) == mol.num_atoms()

    for i in range(0, len(coords)):
        assert coords[i] == mol.atom(i).property("coordinates")

    assert coords[0:3] == [coords[0], coords[1], coords[2]]


def test_atomcharges(ala_mols):
    mols = ala_mols

    mol = mols[0]

    charges = mol.property("charge")

    assert len(charges) == mol.num_atoms()

    for i in range(0, len(charges)):
        assert charges[i] == mol.atom(i).property("charge")

    assert charges[0:3] == [charges[0], charges[1], charges[2]]


def test_atomljs(ala_mols):
    mols = ala_mols

    mol = mols[0]

    ljs = mol.property("LJ")

    assert len(ljs) == mol.num_atoms()

    for i in range(0, len(ljs)):
        assert ljs[i] == mol.atom(i).property("LJ")

    assert ljs[0:3] == [ljs[0], ljs[1], ljs[2]]
