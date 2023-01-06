

def test_trajectory(ala_mols):
    mols = ala_mols
    mol = mols[0]

    assert mol.num_frames() == 0

    assert not mol.has_property("trajectory")

    mol.save_frame()

    assert mol.num_frames() == 1

    assert mol.has_property("trajectory")

    mol.delete_frame(0)

    assert mol.num_frames() == 0

    assert not mol.has_property("trajectory")

    mol.save_frame()

    assert mol.num_frames() == 1

    assert mol.has_property("trajectory")

    mol.save_frame()

    assert mol.num_frames() == 2

    traj = mol.property("trajectory")

    assert traj[0] == traj[1]

    assert traj[0].coordinates() == mol.property("coordinates").to_vector()
