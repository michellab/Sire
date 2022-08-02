
import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    mol = mol.cursor().atoms("element O").apply(
                    lambda a: a.set("special", True)).commit()

    assert len(mol.property("special")) == mol.num_atoms()

    for atom in mol["atom property special"]:
        assert atom.element().num_protons() == 8

    mol = mol.cursor().atoms("element O").delete("special").commit()

    with pytest.raises(KeyError):
        for atom in mol["atom property special"]:
            assert atom.element().num_protons() == 8

    assert len(mol.property("special")) == mol.num_atoms()

    mol = mol.cursor().delete("special").commit()

    with pytest.raises(KeyError):
        p = mol.property("special")


def test_cursor_dict(ala_mols):
    mols = ala_mols
    mol = mols[0]

    cursor = mol.cursor().atoms()

    assert len(cursor) == mol.num_atoms()

    for i, atom in enumerate(mol.atoms()):
        assert cursor[i].name == atom.name().value()
        assert cursor[i].number == atom.number().value()
        assert cursor[i].index == atom.index().value()

        assert cursor[i]["coordinates"] == atom.coordinates()
        assert cursor[i]["charge"] == atom.charge()
        assert cursor[i]["mass"] == atom.mass()

    for atom in cursor.parent().atoms("element O"):
        atom["coordinates"] = (1, 2, 3)

    mol = cursor.parent().commit()

    for atom in mol.atoms("element O"):
        assert atom.x() == 1.0
        assert atom.y() == 2.0
        assert atom.z() == 3.0

    for atom in cursor.parent().atoms("element O"):
        del atom["coordinates"]

    mol = cursor.parent().commit()

    from sire.maths import Vector

    for atom in mol.atoms("element O"):
        assert atom.coordinates() == Vector(0)


def test_cursors(ala_mols):
    mols = ala_mols
    mol = mols[0]

    cursors = mol.cursor().atoms()

    assert len(cursors) == mol.num_atoms()

    assert len(cursors[0:5]) == 5
    assert len(cursors[[0, 2, 4]]) == 3

    assert cursors[-1].name == mol.atom(-1).name().value()

    idxs = range(0,5)
    cs = cursors[0:5]

    for i in range(0,5):
        assert cs[i].id() == cursors[idxs[i]].id()

    idxs = [0, 2, 4]
    cs = cursors[idxs]

    for i in range(0, 3):
        assert cs[i].id() == cursors[idxs[i]].id()
