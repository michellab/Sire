
import pytest


@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


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
