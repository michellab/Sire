
import sire as sr


def test_dcd():
    mols = sr.load(sr.expand(sr.tutorial_url, ["h7n9.pdb", "h7n9.dcd"]))
    mol = mols[0]

    assert mol.num_frames() > 1