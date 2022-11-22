
import pytest


def test_search_property_exists(ala_mols):
    mols = ala_mols

    assert mols["property bond"] == mols.molecules()

    with pytest.raises(KeyError):
        mols["property bnd"]

    mol = mols[0]

    assert mol["property charge"] == mol
    assert mol["atom property charge"] == mol.atoms()

    cursor = mol.cursor()
    cursor["find_me"] = ["complex", "property"]
    mol = cursor.commit()

    assert mol["property find_me"] == mol

    with pytest.raises(KeyError):
        mol["atom property find_me"]

    with pytest.raises(KeyError):
        mols["property find_me"]

    mols.update(mol)

    assert mols["property find_me"] == mol

    with pytest.raises(KeyError):
        mols["atom property find_me"]

    cursor[5]["find_atom"] = True

    with pytest.raises(KeyError):
        mol["property find_atom"]

    mol = cursor.commit()

    assert mol["property find_atom"] == mol
    assert mol["atom property find_atom"] == mol[5]

    assert mol["atom property find_atom == 1"] == mol["atomidx 5"]
    assert mol["atom property find_atom != 1"] == mol["not atomidx 5"]

    with pytest.raises(KeyError):
        mols["property find_atom"]

    mols.update(mol)

    assert mols["property find_atom"] == mol
    assert mols["atom property find_atom"] == mol[5]

    assert mols["property find_me"] == mol

    cursor[4]["find_atom"] = 4

    mol = cursor.commit()

    assert mol["atom property find_atom == 1"] == mol[5]
    assert mol["atom property find_atom == 4"] == mol[4]

    mols.update(mol)

    assert mols["atom property find_atom == 1"] == mol[5]
    assert mols["atom property find_atom == 4"] == mol[4]

    cursor["metadata"] = "something"

    with pytest.raises(KeyError):
        mol["property metadata"]

    mol = cursor.commit()

    assert mol["property metadata"] == mol

    with pytest.raises(KeyError):
        mols["property metadata"]

    mols.update(mol)

    assert mols["property metadata"] == mol

    assert mols["property metadata == something"] == mol

    with pytest.raises(KeyError):
        mols["property metadata == else"]



