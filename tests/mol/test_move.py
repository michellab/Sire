
import pytest


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


def test_cursor_bond_moves(ala_mols):
    mols = ala_mols

    mol = mols[0]

    cursor = mol.cursor()

    lengths0 = cursor.bonds().lengths()

    cursor.bonds().change_length(0.5)

    lengths1 = cursor.bonds().lengths()

    for l0, l1 in zip(lengths0, lengths1):
        assert l0.value() == pytest.approx(l1.value() - 0.5)

    mol = cursor.commit()

    lengths1 = mol.bonds().lengths()

    for l0, l1 in zip(lengths0, lengths1):
        assert l0.value() == pytest.approx(l1.value() - 0.5)

    mol = mols[0]

    cursor = mol.cursor()

    hbonds = cursor.bonds("element H")

    lengths_h0 = hbonds.lengths()
    lengths0 = hbonds.invert().lengths()

    hbonds.change_length(0.5)

    lengths_h1 = hbonds.lengths()
    lengths1 = hbonds.invert().lengths()

    for l0, l1 in zip(lengths0, lengths1):
        assert l0.value() == pytest.approx(l1.value())

    for l0, l1 in zip(lengths_h0, lengths_h1):
        assert l0.value() == pytest.approx(l1.value() - 0.5)

    cursor = mols.cursor()

    lengths = cursor.bonds("element C", "element H").set_lengths(100.0).lengths()

    for l in lengths:
        assert l.value() == pytest.approx(100.0)

    mols = cursor.commit()

    mol = mols[0]

    bonds = mol.bonds("element C", "element H")

    for l in bonds.lengths():
        assert l.value() == pytest.approx(100.0)

    for l in bonds.invert().lengths():
        assert l.value() != pytest.approx(100.0)


def test_cursor_angle_moves(ala_mols):
    mols = ala_mols

    mol = mols[0]

    cursor = mol.cursor()

    cursor.angles()[0].set_size(25)

    mol = cursor.commit()

    angle = mol.angles()[0]

    from sire.units import degrees

    assert angle.size().to(degrees) == pytest.approx(25)

    for s in angle.invert().sizes():
        assert s.to(degrees) != pytest.approx(25)

    cursor = mols.cursor()

    cursor["water"][0:10].angles(
                "element H", "element O", "element H").change_size(10)

    mols = cursor.commit()

    ref_water = mols["water"][10]
    ref_angle = ref_water.angle("element H", "element O", "element H").size()

    for s in mols["water"][0:10].angles(
            "element H", "element O", "element H").sizes():
        assert s.to(degrees) == pytest.approx(ref_angle.to(degrees) + 10.0)


def test_cursor_dihedral_moves(ala_mols):
    mols = ala_mols

    mol = mols[0]

    cursor = mol.cursor()

    # should change all around this dihedral
    cursor.dihedrals()[0].change_size(60)

    mol2 = cursor.commit()

    from sire.units import degrees

    def wrap(ang):
        while ang > 360:
            ang -= 360

        return ang

    for i in range(0,6):
        assert mol2.dihedrals()[i].size().to(degrees) == pytest.approx(
                wrap(mol.dihedrals()[i].size().to(degrees) + 60))

    for i in range(6, len(mol.dihedrals())):
        assert mol2.dihedrals()[i].size().to(degrees) == pytest.approx(
                mol.dihedrals()[i].size().to(degrees))

    mol = mols[0]

    # this should change this dihedral only
    cursor = mol.cursor()

    cursor.dihedrals()[0].change_size(60, move_all=False,
                                      weighting="absolute_mass")

    mol2 = cursor.commit()

    print(mol.dihedrals().sizes())
    print(mol2.dihedrals().sizes())

    for i in [0,3]:
        assert mol2.dihedrals()[i].size().to(degrees) == pytest.approx(
                wrap(mol.dihedrals()[i].size().to(degrees) + 60))

    for i in [1,2] + list(range(4,len(mol.dihedrals()))):
        assert mol2.dihedrals()[i].size().to(degrees) == pytest.approx(
                mol.dihedrals()[i].size().to(degrees))
