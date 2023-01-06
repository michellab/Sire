

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


def test_selector_bonds(ala_mols):
    mols = ala_mols

    mol = mols[0]

    assert len(mol.bonds()) == 21
    assert len(mol["resnum 1"].bonds()) == 5

    for bond in mol["resnum 1"].bonds():
        assert len(bond) == 2
        for atom in bond:
            assert atom.residue().number().value() == 1

    bonds = mol.bonds("element C", "element C")

    assert len(bonds) == 3

    for bond in bonds:
        for atom in bond:
            assert atom.element().num_protons() == 6

    bonds = mol.bonds("element C", "*")

    assert len(bonds) == 19

    for bond in bonds:
        # at least one carbon
        assert bond[0].element().num_protons() == 6 or \
               bond[1].element().num_protons() == 6

    bonds = mol.bonds("element C", "not element C")

    assert len(bonds) == 16

    for bond in bonds:
        # only one carbon
        if bond[0].element().num_protons() == 6:
            assert bond[1].element().num_protons() != 6
        else:
            assert bond[1].element().num_protons() == 6

    bonds = bonds.invert()

    assert len(bonds) == 5

    for bond in bonds:
        for i in range(0, 2):
            assert bond[i] == bond[bond.id()[i]]

        for i, atom in enumerate(bond):
            assert atom == bond[bond.id()[i]]

        # either both carbon or no carbon
        if bond[0].element().num_protons() == 6:
            assert bond[1].element().num_protons() == 6
        else:
            assert bond[1].element().num_protons() != 6


def test_bond_cursor(ala_mols):
    mols = ala_mols
    mol = mols[0]

    c = mol.cursor()

    assert len(c.bonds()) == len(mol.bonds())

    for i, bond in enumerate(c.bonds()):
        bond["count"] = i
        bond["count_string"] = f"string_{i}"

    mol = c.commit()

    for i, bond in enumerate(mol.bonds()):
        assert bond.property("count").value() == i
        assert bond.property("count_string").value() == f"string_{i}"
