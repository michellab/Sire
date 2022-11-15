
import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_search_order(ala_mols):
    mols = ala_mols

    cursor = mols.cursor()

    cursor[0]["is_perturbable"] = True
    cursor[1]["is_perturbable"] = True

    mols = cursor.commit()

    mol0 = mols[0]
    mol1 = mols[1]

    result = mols["property is_perturbable"]

    assert len(result) == 2
    assert result[0] == mol0
    assert result[1] == mol1

    reversed_mols = mols[-1::-1]

    result = reversed_mols["property is_perturbable"]

    assert len(result) == 2
    assert result[0] == mol1
    assert result[1] == mol0

    mols._system.remove(mol0.number())
    mols._molecules = None

    result = mols["property is_perturbable"]

    assert result == mol1

    from sire.legacy.Mol import MGName
    mols._system.add(mol0, MGName("all"))
    mols._molecules = None

    result = mols["property is_perturbable"]

    assert len(result) == 2
    assert result[0] == mol1
    assert result[1] == mol0

    print(mols._system.search("property is_perturable"))

    raise ValueError()

