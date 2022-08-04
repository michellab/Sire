

import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


def test_atomcoords(ala_mols):
    mols = ala_mols

    mol = mols[0]

    coords = mol.property("coordinates")

    assert len(coords) == mol.num_atoms()

    for i in range(0, len(coords)):
        assert coords[i] == mol.atom(i).property("coordinates")

    assert coords[0:3] == [coords[0], coords[1], coords[2]]


def test_atomcharges(ala_mols):
    mols = ala_mols

    mol = mols[0]

    charges = mol.property("charge")

    assert len(charges) == mol.num_atoms()

    for i in range(0, len(charges)):
        assert charges[i] == mol.atom(i).property("charge")

    assert charges[0:3] == [charges[0], charges[1], charges[2]]


def test_atomljs(ala_mols):
    mols = ala_mols

    mol = mols[0]

    ljs = mol.property("LJ")

    assert len(ljs) == mol.num_atoms()

    for i in range(0, len(ljs)):
        assert ljs[i] == mol.atom(i).property("LJ")

    assert ljs[0:3] == [ljs[0], ljs[1], ljs[2]]


def test_atomelements(ala_mols):
    mols = ala_mols

    for e in mols["element C"].property("element"):
        assert e.num_protons() == 6

    for e in mols["element H"].apply("element"):
        assert e.num_protons() == 1

    for e in mols["element O"].apply("property", "element"):
        assert e.num_protons() == 8

    for e in mols["element N"].apply(lambda atom: atom.element()):
        assert e.num_protons() == 7

    for e in mols["element C"].apply(lambda atom, key : atom.property(key), "element"):
        assert e.num_protons() == 6

    with pytest.raises(AttributeError):
        mols["element C"].element()


def test_convenience_atom_funcs(ala_mols):
    mols = ala_mols
    mol = mols[0]

    from sire.units import g_per_mol, mod_electron

    mass = 0 * g_per_mol

    for i in range(0, 5):
        mass += mol[i].mass()

    assert mass.approx_equal(mol[0:5].mass())

    charge = 0 * mod_electron

    for atom in mol["resnum 1"].atoms():
        charge += atom.charge()

    assert charge.approx_equal(mol["resnum 1"].charge())

    assert mol.coords() == mol.evaluate().center_of_mass()


def test_atompropprop(ala_mols):
    mols = ala_mols
    mol = mols[0]

    import sire as sr

    cursor = mol.cursor()

    for atom in cursor.atoms():
        atom["charge2"] = 5 * sr.units.mod_electron
        atom["mass2"] = 3 * sr.units.g_per_mol
        atom["coords2"] = sr.maths.Vector(1, 2, 3)

    mol = cursor.commit()

    # all of the above properties should *not* be AtomPropertyProperty...
    assert type(mol.property("charge2")) == type(mol.property("charge"))
    assert type(mol.property("mass2")) == type(mol.property("mass"))
    assert type(mol.property("coords2")) == type(mol.property("coordinates"))
