

def test_selector_dihedrals(ala_mols):
    mols = ala_mols

    mol = mols[0]

    assert len(mol.dihedrals()) == 41
    assert len(mol["resnum 1"].dihedrals()) == 3

    for dihedral in mol["resnum 1"].dihedrals():
        assert len(dihedral) == 4
        for atom in dihedral:
            assert atom.residue().number().value() == 1

    dihedrals = mol.dihedrals("not element C", "element C",
                              "element C", "not element C")

    assert len(dihedrals) == 16

    for dihedral in dihedrals:
        assert dihedral[0].element().num_protons() != 6
        assert dihedral[1].element().num_protons() == 6
        assert dihedral[2].element().num_protons() == 6
        assert dihedral[3].element().num_protons() != 6

    dihedrals = mol.dihedrals("*", "element C", "element C", "*")

    assert len(dihedrals) == 21

    for dihedral in dihedrals:
        assert dihedral[1].element().num_protons() == 6
        assert dihedral[2].element().num_protons() == 6

    dihedrals = dihedrals.invert()

    assert len(dihedrals) == 20

    for dihedral in dihedrals:
        for i in range(0, 4):
            assert dihedral[i] == dihedral[dihedral.id()[i]]

        for i, atom in enumerate(dihedral):
            assert atom == dihedral[dihedral.id()[i]]

        # are allowed to have X-C-X-X or X-X-C-X

        if dihedral[1].element().num_protons() == 6:
            assert dihedral[2].element().num_protons() != 6


def test_dihedral_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    c = mol.cursor()

    assert len(c.dihedrals()) == len(mol.dihedrals())

    for i, dihedral in enumerate(c.dihedrals()):
        dihedral["count"] = i
        dihedral["count_string"] = f"string_{i}"

    mol = c.commit()

    for i, dihedral in enumerate(mol.dihedrals()):
        assert dihedral.property("count").value() == i
        assert dihedral.property("count_string").value() == f"string_{i}"
