
import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")

def test_string_save(ala_mols):
    mols = ala_mols

    from sire import save_to_string

    #s = save_to_string(mols, format="pdb")

    # check that we can load this from a string
    from sire.legacy.IO import PDB2

    #mols2 = PDB2(s).to_system()

    #assert mols.num_atoms() == mols2.num_atoms()
    #assert mols.num_residues() == mols2.num_residues()

    #for i in range(0, len(mols[0])):
    #    assert mols[0][i].name() == mols2[0][i].name()

    # check that we can save a sub-view to a string
    s = save_to_string(mols["element O"], format="pdb")
