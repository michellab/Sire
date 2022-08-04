
import pytest

@pytest.fixture(scope="session")
def chol_mols():
    import sire as sr
    return sr.load_test_files("cholesterol.sdf")


def test_bond_props(chol_mols):
    mols = chol_mols
    mol = mols[0]

    for bond in mol.bonds():
        assert bond.has_property("type")
        assert bond.has_property("sdf_fields")
        assert bond.has_property("stereoscopy")

    for bond in mol.cursor().bonds():
        assert bond["type"] == bond.view().property("type")
        assert bond["sdf_fields"] == bond.view().property("sdf_fields")
        assert bond["stereoscopy"] == bond.view().property("stereoscopy")
