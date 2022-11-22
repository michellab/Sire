
import sire as sr
import pytest

@pytest.fixture(scope="session")
def ala_mols():
    return sr.load_test_files("ala.top", "ala.crd")


@pytest.fixture(scope="session")
def h7n9_mols():
    return sr.load_test_files("h7n9.pdb", "h7n9.dcd")


@pytest.fixture(scope="session")
def chol_mols():
    return sr.load_test_files("cholesterol.sdf")


@pytest.fixture(scope="session")
def ala_traj():
    return sr.load_test_files("ala.top", "ala.traj")


@pytest.fixture(scope="session")
def p38_mols():
    return sr.load_test_files("p38.pdb")


@pytest.fixture(scope="session")
def alanin_mols():
    return sr.load_test_files("alanin.psf")
