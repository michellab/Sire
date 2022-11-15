
import pytest

from sire.legacy.IO import PDB2, Mol2, AmberPrm, GroTop, SDF, CharmmPSF

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")

@pytest.mark.parametrize("format, Loader",
                          [
                              ("pdb", PDB2),
                              #("prmtop", AmberPrm),  This fails as AmberPrm expects full molecules
                              #("mol2", Mol2),        This fails as Mol2 expects full molecules
                              #("top", GroTop),       This fails as GroTop expects full molecules
                              #("sdf", SDF),          This fails as SDF expects full molecules
                              #("psf", CharmmPSF),    This fails as CharmmPSF expects full molecules
                          ])
def test_string_save(ala_mols, format, Loader):
    mols = ala_mols

    from sire import save_to_string

    s = save_to_string(mols, format=format)

    # check that we can load this from a string
    mols2 = Loader(s).to_system()

    assert mols.num_atoms() == mols2.num_atoms()

    if format == "sdf":
        # SDF doesn't save residue information, so we lose
        # that the first molecule has 3 residues
        assert mols.num_residues()-2 == mols2.num_residues()
    else:
        assert mols.num_residues() == mols2.num_residues()

    if format == "sdf":
        # SDF doesn't save the atom names...
        for i in range(0, len(mols[0])):
            assert mols[0][i].element() == mols2[0][i].element()
    else:
        for i in range(0, len(mols[0])):
            assert mols[0][i].name() == mols2[0][i].name()

    # check that we can save a sub-view to a string
    subview = mols["element O"]
    s = save_to_string(mols["element O"], format=format)

    mols2 = Loader(s).to_system()

    assert subview.num_atoms() == mols2.num_atoms()
    assert subview.num_residues() == mols2.num_residues()

    if format == "sdf":
        for i in range(0, len(mols[0])):
            assert mols[0][i].element() == mols2[0][i].element()
    else:
        for i in range(0, len(subview[0])):
            assert subview[0][i].name() == mols2[0][i].name()
