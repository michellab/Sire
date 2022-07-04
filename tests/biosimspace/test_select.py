
import sire as sr
sr.use_mixed_api()

import sire.legacy.Mol as _Mol

import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_selector_wrap(ala_mols):
    mols = ala_mols

    # These will fail if importing Sire under a mixed
    # api and then importing sire.legacy modules causes
    # an issue
    print(mols["atoms within 2 of element C"])

    res = mols["resnum 1"]
