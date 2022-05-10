
import sire as sr

def _assert_same_bonds(b0, b1):
    assert len(b0) == len(b1)

    for bond0 in b0:
        same = False

        for bond1 in b1:
            if bond0 == bond1:
                same = True
                break

        if not same:
            print(f"Not the same!\n{b0}\n{b1}")
            assert False


def test_index_bonds():
    mols = sr.load_test_files("ala.top", "ala.crd")

    mol = mols[0]

    assert len(mol["bonds in resnum 1"]) == 5
    assert len(mol["bonds involving resnum 1"]) == 6

    _assert_same_bonds(mol["bonds in resnum 1"], mol["resnum 1"].bonds())
    _assert_same_bonds(mol["bonds involving resnum 1"], mol.bonds("resnum 1"))
    _assert_same_bonds(mol["bonds in resnum 2"], mol["bonds within resnum 2"])


    cx = mol["bonds to element C"]
    ccx = mol["bonds involving element C"]
    cc = mol["bonds from element C to element C"]

    assert len(cc) + len(cx) == len(ccx)

    assert len(cc) == 3
    assert len(cx) == 16
    assert len(ccx) == 19

    _assert_same_bonds(ccx["bonds from element C to element C"], cc)


def test_index_mols_bonds():
    mols = sr.load_test_files("ala.top", "ala.crd")

    bnds = mols["bonds from element O to element H"]

    assert len(bnds) == 2 * (mols.num_molecules() - 1)
    assert len(bnds) == 1260  # just in case something else failed!

    assert bnds.num_molecules() == mols.num_molecules() - 1
    assert bnds.num_atoms() == bnds.num_molecules() * 3


if __name__ == "__main__":
    test_index_bonds()
    test_index_mols_bonds()


