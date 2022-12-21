
import pytest


def assert_approx_equal(nrg0, nrg1):
    assert nrg0.value() == pytest.approx(nrg1.value())

    for component in nrg0.components():
        try:
            assert nrg0[component].value() == pytest.approx(nrg1[component].value())
        except Exception as e:
            print(f"FAILED {component}")
            print(nrg0.components())
            print(nrg1.components())
            raise e


def test_trajectory_energies(ala_traj):
    mols = ala_traj
    mol = mols[0]

    from sire.units import angstrom
    from sire.base import create_map

    map = create_map({"cutoff": 10*angstrom})

    def assert_same(view):
        from sire import colname

        nrgs0 = [0, 0, 0]

        view.load_frame(0)
        nrgs0[0] = view.energies(map=map)
        view.load_frame(1)
        nrgs0[1] = view.energies(map=map)
        view.load_frame(2)
        nrgs0[2] = view.energies(map=map)

        nrgs1 = view.trajectory(map=map)[0:3].energies(to_pandas=False)

        for idx, v in enumerate(view):
            for i in range(0, 3):
                assert nrgs0[i][idx].value() == pytest.approx(nrgs1[colname(v, "total")][i])

    assert_same(mol.residues())

    def assert_same_pair(view, other):
        from sire import colname

        nrgs0 = [0, 0, 0]

        view.load_frame(0)
        other.load_frame(0)
        nrgs0[0] = view.energies(other, map=map)
        view.load_frame(1)
        other.load_frame(1)
        nrgs0[1] = view.energies(other, map=map)
        view.load_frame(2)
        other.load_frame(2)
        nrgs0[2] = view.energies(other, map=map)

        nrgs1 = view.trajectory()[0:3].energies(other, to_pandas=False, map=map)

        for idx, v in enumerate(view):
            for i in range(0, 3):
                assert nrgs0[i][idx].value() == pytest.approx(nrgs1[colname(v, "total")][i], 1e-5)

    assert_same_pair(mol.residues()[0:2], mol.residues()[-1])
    assert_same_pair(mols["water"][0:5], mol)
    assert_same_pair(mol.atoms(), mols["water and molecule within 5 of molidx 0"])


def test_trajectory_energy(ala_traj):
    mols = ala_traj
    mol = mols[0]

    from sire.units import angstrom
    from sire.base import create_map

    map = create_map({"cutoff": 10*angstrom})

    def assert_same(view):
        view.load_frame(0)
        nrg00 = view.energy(map=map)
        view.load_frame(1)
        nrg01 = view.energy(map=map)
        view.load_frame(2)
        nrg02 = view.energy(map=map)

        nrgs = view.trajectory()[0:3].energy(to_pandas=False, map=map)

        def check(n0, n1, i):
            if n0.value() != pytest.approx(n1["total"][i]):
                print(i, n0.value(), n1["total"][i])
                for key in n0.components().keys():
                    print(key, n0[key], n1[key])

            assert n0.value() == pytest.approx(n1["total"][i])

        check(nrg00, nrgs, 0)
        check(nrg01, nrgs, 1)
        check(nrg02, nrgs, 2)

    assert_same(mol)
    assert_same(mol.bonds()[0])
    assert_same(mol.atoms("element C"))
    assert_same(mols)

    def assert_same_pair(view, other):
        view.load_frame(0)
        other.load_frame(0)
        nrg00 = view.energy(other, map=map)
        view.load_frame(1)
        other.load_frame(1)
        nrg01 = view.energy(other, map=map)
        view.load_frame(2)
        other.load_frame(2)
        nrg02 = view.energy(other, map=map)

        nrgs = view.trajectory()[0:3].energy(other, to_pandas=False,
                                             map=map)

        def check(n0, n1, i):
            if n0.value() != pytest.approx(n1["total"][i]):
                print(i, n0.value(), n1["total"][i])
                for key in n0.components().keys():
                    print(key, n0[key], n1[key])

            assert n0.value() == pytest.approx(n1["total"][i])

        check(nrg00, nrgs, 0)
        check(nrg01, nrgs, 1)
        check(nrg02, nrgs, 2)

    assert_same_pair(mol["element C"], mol["not element C"])
    assert_same_pair(mol, mols["water"])
    assert_same_pair(mol.residues()[0:2], mol.residues()[-1])


def test_energy(ala_mols):
    mols = ala_mols
    mol = mols[0]

    total = mol.energy()

    assert_approx_equal(total, mol.atoms().energy())
    assert_approx_equal(total, mol.residues().energy())
    assert_approx_equal(total, mols[0].energy())
    assert_approx_equal(total, mols["molidx 0"].energy())

    assert total["bond"].value() == pytest.approx(mol.bonds().energy().value())
    assert total["angle"].value() == pytest.approx(mol.angles().energy().value())
    assert total["dihedral"].value() == pytest.approx(mol.dihedrals().energy().value())
    assert total["improper"].value() == pytest.approx(mol.impropers().energy().value(), rel=1e-5)

    total0 = mol["element C"].energy()
    total1 = mol["not element C"].energy()
    total2 = mol["element C"].energy(mol["not element C"])

    assert_approx_equal(total, total0 + total1 + total2)

    total2 = mol["not element C"].energy(mol["element C"])

    assert_approx_equal(total, total0 + total1 + total2)

    assert total0["bond"].value() == pytest.approx(mol.bonds("element C", "element C").energy().value())
    assert total1["bond"].value() == pytest.approx(mol.bonds("not element C", "not element C").energy().value())
    assert total2["bond"].value() == pytest.approx(mol.bonds("element C", "not element C").energy().value())

    total = mol.residues()[0:2].energy()

    total0 = mol.residues()[0].energy()
    total1 = mol.residues()[1].energy()
    total2 = mol.residues()[0].energy(mol.residues()[1])

    assert_approx_equal(total, total0 + total1 + total2)

    total2 = mol.residues()[1].energy(mol.residues()[0])

    assert_approx_equal(total, total0 + total1 + total2)

    total = mol["atomidx < 6"].energy()

    total0 = mol["atomidx < 3"].energy()
    total1 = mol["atomidx >=3 and atomidx < 6"].energy()
    total2 = mol["atomidx < 3"].energy(mol["atomidx >=3 and atomidx < 6"])

    assert_approx_equal(total, total0 + total1 + total2)

    total2 = mol["atomidx >=3 and atomidx < 6"].energy(mol["atomidx < 3"])

    assert_approx_equal(total, total0 + total1 + total2)

    total = mols["water"].energy()

    assert total["angle"] == pytest.approx(mols["water"].angles().energy())

    total = mols.energy()

    total0 = mols[0].energy()
    total1 = mols["water"].energy()
    total2 = mols[0].energy(mols["water"])

    assert_approx_equal(total, total0 + total1 + total2)

    total2 = mols["water"].energy(mols[0])

    assert_approx_equal(total, total0 + total1 + total2)

    total = mols["element O"].energy()

    total0 = mols["water and element O"].energy()
    total1 = mols["(not water) and element O"].energy()
    total2 = mols["water and element O"].energy(mols["(not water) and element O"])

    assert_approx_equal(total, total0 + total1 + total2)

    total2 = mols["(not water) and element O"].energy(mols["water and element O"])


def test_neura_energy(neura_mols):
    mols = neura_mols

    from sire.units import angstrom

    # this is an integration test, plus test that 'map' is
    # being interpreted and passed correctly
    components = mols[0:5].energy(map={"cutoff": 15*angstrom}).components()

    # these values have been pre-calculated. The test checks
    # if anything has changed the energies
    assert components["coulomb"].value() == pytest.approx(-587.683)
    assert components["LJ"].value() == pytest.approx(-17.15338212)

    components = mols[0:5].energy(map={"cutoff": 5*angstrom}).components()

    # these values have been pre-calculated. The test checks
    # if anything has changed the energies
    assert components["coulomb"].value() == pytest.approx(-223.95138492)
    assert components["LJ"].value() == pytest.approx(2.09059, 1e-5)
