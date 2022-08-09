import sire as sr
import pytest


@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_selector_impropers(ala_mols):
    mols = ala_mols

    mol = mols[0]

    assert len(mol.impropers()) == 4
    assert len(mol["resnum 1"].impropers()) == 0

    assert len(mol.impropers("resnum 1", "resnum 2")) == 2
    assert len(mol.impropers("resnum 2", "resnum 3")) == 2

    impropers = mol.impropers("element C", "*", "element C", "*")

    assert len(impropers) == 2

    for improper in impropers:
        assert improper[0].element().num_protons() == 6
        assert improper[2].element().num_protons() == 6

    impropers = impropers.invert()

    assert len(impropers) == 2

    for improper in impropers:
        for i in range(0, 4):
            assert improper[i] == improper[improper.id()[i]]

        for i, atom in enumerate(improper):
            assert atom == improper[improper.id()[i]]

        if improper[0].element().num_protons() == 6:
            assert improper[2].element().num_protons() != 6


def test_improper_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    c = mol.cursor()

    assert len(c.impropers()) == len(mol.impropers())

    for i, improper in enumerate(c.impropers()):
        improper["count"] = i
        improper["count_string"] = f"string_{i}"

    mol = c.commit()

    for i, improper in enumerate(mol.impropers()):
        assert improper.property("count").value() == i
        assert improper.property("count_string").value() == f"string_{i}"
