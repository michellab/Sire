
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


def test_cursor_assign_from_container(ala_mols):
    mols = ala_mols
    mol = mols[0]

    cursor = mol.cursor()

    cursor["from_dict"] = {"cat": "meow", "dog": "woof", "fish": 3}

    mol = cursor.commit()

    assert mol.property("from_dict")["cat"] == "meow"
    assert mol.property("from_dict")["dog"] == "woof"
    assert mol.property("from_dict")["fish"] == 3

    cursor["from_list"] = ["cat", "dog", "horse"]

    mol = cursor.commit()

    assert mol.property("from_list")[0] == "cat"
    assert mol.property("from_list")[1] == "dog"
    assert mol.property("from_list")[2] == "horse"

    cursor["from_list"] = ["cat", 5]

    mol = cursor.commit()

    assert mol.property("from_list")[0] == "cat"
    assert mol.property("from_list")[1] == 5

    cursor = cursor[0]

    cursor["atom_test"] = ["dog", 3.1, 5]

    mol = cursor.parent().commit()

    assert mol[0].property("atom_test")[0] == "dog"
    assert mol[0].property("atom_test")[1] == 3.1
    assert mol[0].property("atom_test")[2] == 5

    cursor = cursor.parent()[1]

    cursor["atom_test"] = {"cat": "meow", "dog": "woof"}

    mol = cursor.parent().commit()

    assert mol[0].property("atom_test")[0] == "dog"
    assert mol[0].property("atom_test")[1] == 3.1
    assert mol[0].property("atom_test")[2] == 5

    assert mol[1].property("atom_test")["cat"] == "meow"
    assert mol[1].property("atom_test")["dog"] == "woof"


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

    from sire.units import angstrom

    for atom in mol.atoms("element O"):
        assert atom.x() == 1.0 * angstrom
        assert atom.y() == 2.0 * angstrom
        assert atom.z() == 3.0 * angstrom

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


def test_cursor_renaming(ala_mols):
    mols = ala_mols
    mol = mols[0]

    cursor = mol.cursor()

    for atom in cursor.atoms():
        atom.name = f"{atom.index}"
        atom.number = atom.index

    mol = cursor.commit()

    for atom in mol.atoms():
        assert atom.name().value() == f"{atom.index().value()}"
        assert atom.number().value() == atom.index().value()

    mol = mols[0]

    for atom in mol.atoms():
        assert atom.name().value() != f"{atom.index().value()}"
        assert atom.number().value() != atom.index().value()

    mol = mol.cursor().atoms().apply(
        lambda atom: atom.set_name(f"{atom.index}")
    ).apply(
        lambda atom: atom.set_number(atom.index)
    ).commit()

    for atom in mol.atoms():
        assert atom.name().value() == f"{atom.index().value()}"
        assert atom.number().value() == atom.index().value()


def test_inverse_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    c0 = mol.cursor().atoms("element C")
    c1 = mol.atoms("element C").cursor()

    assert len(c0) == len(c1)

    for i in range(0, len(c0)):
        assert c0[i].id() == c1[i].id()

    assert c0.commit().number() == c1.commit().number()

    c0 = mol.cursor().bonds()
    c1 = mol.bonds().cursor()

    assert len(c0) == len(c1)

    for i in range(0, len(c0)):
        assert c0[i].id() == c1[i].id()

    assert c0.commit().number() == c1.commit().number()


def test_cursor_translate(ala_mols):
    mols = ala_mols
    mol = mols[0]

    cursor = mol.cursor()

    coords = cursor["coordinates"]

    cursor.translate(1, 0, 0)

    for i in range(0, len(coords)):
        c0 = coords[i]
        c1 = cursor["coordinates"][i] - (1,0,0)
        assert c0[0].value() == pytest.approx(c1[0].value())
        assert c0[1] == c1[1]
        assert c0[2] == c1[2]

    cursor.translate((2,3,4))

    for i in range(0, len(coords)):
        c0 = coords[i]
        c1 = cursor["coordinates"][i] - (3,3,4)
        assert c0[0].value() == pytest.approx(c1[0].value())
        assert c0[1].value() == pytest.approx(c1[1].value())
        assert c0[2].value() == pytest.approx(c1[2].value())

    cursor = mol["element C"].cursor()

    coords = mol.property("coordinates")

    cursor.translate((2,3,4))

    mol = cursor.commit()

    from sire.mol import Element
    carbon = Element("Carbon")

    for i in range(0, len(coords)):
        atom = mol[i]
        c0 = coords[i]

        if atom.element() == carbon:
            c1 = atom.coordinates() - (2,3,4)
        else:
            c1 = atom.coordinates()

        assert c0[0].value() == pytest.approx(c1[0].value())
        assert c0[1].value() == pytest.approx(c1[1].value())
        assert c0[2].value() == pytest.approx(c1[2].value())

    cursor = mols.cursor()

    from sire.maths import Vector
    cursor.translate(Vector(2, 3, 4), map={"coordinates":"coordinates"})

    mol100 = mols[100]
    c100 = cursor[100]

    coords = mol100.property("coordinates")

    for i in range(0, len(coords)):
        c0 = coords[i]
        c1 = c100["coordinates"][i] - (2, 3, 4)

        assert c0[0].value() == pytest.approx(c1[0].value())
        assert c0[1].value() == pytest.approx(c1[1].value())
        assert c0[2].value() == pytest.approx(c1[2].value())

    mols = cursor.commit()

    coords1 = mols[100].property("coordinates")

    for i in range(0, len(coords)):
        c0 = coords[i]
        c1 = coords1[i] - (2, 3, 4)

        assert c0[0].value() == pytest.approx(c1[0].value())
        assert c0[1].value() == pytest.approx(c1[1].value())
        assert c0[2].value() == pytest.approx(c1[2].value())


def test_cursor_indexing(ala_mols):
    mols = ala_mols

    cursor = mols.cursor()

    assert len(cursor) == len(mols)

    c = cursor[0]

    assert c.number == mols[0].number().value()
    assert c.view().selected_all()

    assert c.is_same_editor(c.atoms())
    assert c.is_same_editor(c.atoms()[0])
    assert c.atoms()[0].is_same_editor(c)

    mol = mols[0]

    assert len(c.atoms()) == len(mol.atoms())
    assert len(c.atoms()[0:8:2]) == len(mol.atoms()[0:8:2])

    assert c.atoms()[0:8:2].view() == mol.atoms()[0:8:2]

    assert c.atoms("element C").view() == mol["element C"]
    assert cursor["element C"].view() == mols["element C"]

    assert cursor["element C"].commit() == mols["element C"]

    assert cursor["element C"][0].is_same_editor(c)

    assert cursor["element C"].molecule().view().selected_all()
    assert cursor["element C"].molecule().is_same_editor(c)

    assert c["element C"]["element"][0] == mol["element C"][0].property("element")

    with pytest.raises(Exception):
        c("element")

    with pytest.raises(KeyError):
        c["elemen"]