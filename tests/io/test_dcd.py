

def test_dcd(h7n9_mols):
    mols = h7n9_mols
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


def test_dcd_trajectory(h7n9_mols):
    mols = h7n9_mols

    i = 0

    for t in mols.trajectory()[0::100]:
        i += 1

    assert i == (mols.num_frames() // 100) + 1

    bonds = mols.bonds()[0:5]

    for t in bonds.trajectory()[0::100]:
        assert t.count() == bonds.count()
        assert len(t.measures()) == len(bonds.measures())

    t = bonds.trajectory()[100].__next__()

    bonds.load_frame(100)

    assert bonds.measures() == t.measures()
