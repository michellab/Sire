import sire as sr
import pytest


@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_selector_angles(ala_mols):
    mols = ala_mols

    mol = mols[0]

    assert len(mol.angles()) == 36
    assert len(mol["resnum 1"].angles()) == 7

    for angle in mol["resnum 1"].angles():
        assert len(angle) == 3
        for atom in angle:
            assert atom.residue().number().value() == 1

    angles = mol.angles("element C", "element C", "element C")

    assert len(angles) == 1

    for atom in angles[0]:
        assert atom.element().num_protons() == 6

    angles = mol.angles("*", "element C", "*")

    assert len(angles) == 30

    for angle in angles:
        assert angle[1].element().num_protons() == 6

    angles = angles.invert()

    assert len(angles) == 6

    for angle in angles:
        for i in range(0, 3):
            assert angle[i] == angle[angle.id()[i]]

        for i, atom in enumerate(angle):
            assert atom == angle[angle.id()[i]]

        assert angle[1].element().num_protons() != 6


def test_angle_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    c = mol.cursor()

    assert len(c.angles()) == len(mol.angles())

    for i, angle in enumerate(c.angles()):
        angle["count"] = i
        angle["count_string"] = f"string_{i}"

    mol = c.commit()

    for i, angle in enumerate(mol.angles()):
        assert angle.property("count").value() == i
        assert angle.property("count_string").value() == f"string_{i}"


def test_multi_angles(ala_mols):
    mols = ala_mols

    from sire.mm import SelectorMAngle

    mol = mols[4]

    angs = SelectorMAngle(mol.atoms("element H"), mol.atoms("element O"), mol.atoms("element H"))

    for ang in angs:
        assert ang[0].element().num_protons() == 1
        assert ang[1].element().num_protons() == 8
        assert ang[2].element().num_protons() == 1

    angs = mols.angles("element H", "element O", "element H")

    for ang in angs:
        assert ang[0].element().num_protons() == 1
        assert ang[1].element().num_protons() == 8
        assert ang[2].element().num_protons() == 1
