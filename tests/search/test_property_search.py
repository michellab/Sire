
import pytest


def test_search_property_map(ala_mols):
    """Tests that searches that involve properties can be mapped to
       different properties by passing in a property map
    """
    mols = ala_mols.clone()

    mol = mols[0]

    cursor = mol.cursor()

    c = cursor["connectivity"]
    c = c.edit().disconnect_all().commit()
    cursor["connectivity2"] = c

    cursor["charge2"] = cursor["charge"]
    cursor["mass2"] = cursor["mass"]

    for i in range(0, len(cursor)):
        cursor[i]["charge2"] *= 1000
        cursor[i]["mass2"] *= 1000

    mol = cursor.commit()

    with pytest.raises(KeyError):
        mol["atom charge > 10"]

    assert mol[("atom charge > 10", {"charge":"charge2"})] == mol["atom charge > 0.01"]
    assert mol[("atom charge < 10", {"charge":"charge2"})] == mol["atom charge < 0.01"]

    assert mol[("atom mass > 10", {"mass":"mass2"})] == mol["atom mass > 0.01"]

    with pytest.raises(KeyError):
        mol[("bonds to element C", {"connectivity":"connectivity2"})]

    assert len(mol["bonds to element C"]) > 0

    with pytest.raises(KeyError):
        mols["property kitten"]

    assert mols[("property kitten", {"kitten":"charge"})] == mols.molecules()

    cursor[5]["find_atom"] = True

    mol = cursor.commit()
    mols.update(mol)

    assert mols["atom property find_atom"] == mol[5]

    with pytest.raises(KeyError):
        mols["atom property find_mapped"]

    assert mols[("atom property find_mapped", {"find_mapped": "find_atom"})] == mol[5]


def test_search_property_exists(ala_mols):
    mols = ala_mols.clone()

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



