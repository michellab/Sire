
import sire as sr


def test_dcd():
    mols = sr.load(sr.expand(sr.tutorial_url, ["h7n9.pdb", "h7n9.dcd"]))
    mol = mols[0]

    assert mol.num_frames() == 501
    assert mols.molecules().num_frames() == 501
    assert mols.atoms().num_frames() == 501

    coords = mol.coordinates()

    mol.load_frame(0)

    assert coords == mol.coordinates()

    mol.load_frame(1)

    assert coords != mol.coordinates()

    coords01 = mol.coordinates()

    mol.load_frame(0)

    assert coords == mol.coordinates()

    mol.load_frame(1)

    assert coords01 == mol.coordinates()

    mol.load_frame(500)

    coords500 = mol.coordinates()

    mol.load_frame(1)

    assert coords01 == mol.coordinates()

    mol.load_frame(500)

    assert coords500 == mol.coordinates()

    molecules = mols.molecules()

    molecules.load_frame(500)

    assert molecules[0].coordinates() == coords500
