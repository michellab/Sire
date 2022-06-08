
import pytest

@pytest.fixture(scope="session")
def ala_mols():
    import sire as sr
    return sr.load_test_files("ala.top", "ala.crd")


@pytest.fixture(scope="session")
def p38_mols():
    import sire as sr
    return sr.load_test_files("p38.pdb")


@pytest.fixture(scope="session")
def alanin_mols():
    import sire as sr
    return sr.load_test_files("alanin.psf")


def test_distance_searching(ala_mols):
    mols = ala_mols

    print(mols["atoms within 2.0 angstrom of element C"])



def test_bond_property_searching(ala_mols):
    mols = ala_mols.clone()

    cursor = mols[0].cursor()

    b = cursor.bonds()[0]

    b["is_perturbable"] = True
    b["test property"] = "cat goes meow"
    b["number"] = 312.1

    mols.update(cursor.commit())

    bond = mols[0].bond(b.id())

    assert mols["bond property is_perturbable"] == bond
    assert mols["bond property 'test property' == 'cat goes meow'"] == bond
    assert mols["bond property number > 312"] == bond


def test_res_property_searching(ala_mols):
    mols = ala_mols.clone()

    res = mols[0]["residx 0"]

    mols.update(res.edit()
                    .set_property("is_perturbable", True)
                    .set_property("test property", "cat goes meow")
                    .set_property("number", -3.5).molecule().commit()
               )

    assert mols["residue property is_perturbable"] == mols[0]["residx 0"]
    assert mols["residue property 'test property' == 'cat goes meow'"] == mols[0]["residx 0"]
    assert mols["residue property number < -3.4"] == mols[0]["residx 0"]


def test_atom_property_searching(ala_mols):
    mols = ala_mols.clone()

    mols.update(mols[0][0].edit().set_property("is_perturbable", True).molecule().commit())
    mols.update(mols[0][1].edit().set_property("test property", "cat goes meow").molecule().commit())
    mols.update(mols[0][2].edit().set_property("is_perturbable", True).molecule().commit())
    mols.update(mols[0][3].edit().set_property("number", 5.4).molecule().commit())

    assert mols["atom property is_perturbable"] == mols[0][ [0,2] ]
    assert mols["atom property 'test property' == 'cat goes meow'"] == mols[0][1]
    assert mols["atom property number > 5.3"] == mols[0][3]
    assert len(mols["atom property number < 5.4"]) == mols[0].num_atoms() - 1



def test_property_searching(ala_mols):
    mols = ala_mols.clone()

    mols.update(mols[0].edit()
                       .set_property("test property", "cat goes meow")
                       .set_property("is_perturbable", True)
                       .set_property("false_property", False)
                       .set_property("number", 5.4)
                       .set_property("val", -42)
                       .commit())
   
    assert mols["property 'test property' == 'cat goes meow'"] == mols[0]

    assert mols["property is_perturbable"] == mols[0]
    assert mols["property is_perturbable == true"] == mols[0]
    assert mols["property false_property == off"] == mols[0]
    assert mols["property number == 5.4"] == mols[0]
    assert mols["property number =~ 5.4"] == mols[0]
    assert mols["property number > 3"] == mols[0]
    assert mols["property number < 6"] == mols[0]
    assert mols["property val == -42"] == mols[0]
    assert mols["property val < 0"] == mols[0]
    assert mols["property val >= -42"] == mols[0]

    with pytest.raises(KeyError):
        mols["property number > 5.4"]

    with pytest.raises(KeyError):
        mols["property false_property == True"]

    with pytest.raises(KeyError):
        mols["property val > 0"]


def test_basic_indexing(ala_mols):
    mols = ala_mols
    import sire as sr

    assert len(mols["atomname CH3"]) == 2
    assert len(mols[sr.atomid("CH3")]) == 2

    with pytest.raises(KeyError):
        mols[sr.atomid("ch3")]

    assert len(mols[sr.atomid("ch3", case_sensitive=False)]) == 2

    with pytest.raises(KeyError):
        mols["atomname ch3"]

    assert len(mols["atomname /ch3/i"]) == 2
    assert len(mols["atomname /CH*/"]) == 2

    assert len(mols["atomidx 0"]) == mols.num_molecules()
    assert len(mols["atomidx 0,1"]) == 2 * mols.num_molecules()
    assert len(mols["atomidx 0:3"]) == 3 * mols.num_molecules()
    assert len(mols["atomidx 0:4"]) == 4 + (3 * (mols.num_molecules()-1))

    assert len(mols["atomnum 1"]) == 1
    assert len(mols["atomnum 1,2"]) == 2
    
    assert len(mols["atoms in *"]) == mols.num_atoms()
    assert len(mols["residues in *"]) == mols.num_residues()


def test_search_indexing(ala_mols):
    mols = ala_mols

    assert mols["{element C}[0]"] == mols["element C"][0]
    assert mols["{element C}[-1]"] == mols["element C"][-1]
    assert mols["{element C}[0:2]"] == mols["element C"][0:2]
    assert mols["{element C}[2:0:-1]"] == mols["element C"][2:0:-1]


def test_logical_indexing(ala_mols):
    mols = ala_mols

    assert len(mols["element C or element O"]) == len(mols["element C"]) + len(mols["element O"])

    assert len(mols["resname ALA or resname ACE"]) == len(mols["resname ALA,ACE"])

    assert len(mols["not resname ALA"]) == mols.num_residues() - 1

    assert len(mols["resname ALA and element C"]) == len(mols["resname ALA"]["element C"])

    assert mols["resname ALA and element C"] == mols["element C and resname ALA"]

    assert len(mols["not molidx 0"]) == mols.num_molecules() - 1

    assert len(mols["atomname HH31 or atomname HH32 or atomname HH33"]) == len(mols["atomname HH31, HH32, HH33"])

    assert len(mols["(element C or element O) and (element O or element H)"]) == len(mols["element O"])

    assert len(mols["not element C"]) == len(mols["element H or element O or element N"])

    assert len(mols["(not element C) and (not element O) and (not element H)"]) == len(mols["element N"])

def test_complex_indexing(p38_mols):
    import sire as sr

    mols = p38_mols
    mol = mols[0]

    assert mols[f"molname {mol.name().value()}"] == mol
    assert mols["molidx 0"] == mol
    assert mols[f"molnum {mol.number().value()}"] == mol

    assert mol.atom(0).index() == sr.mol.AtomIdx(0)
    assert mol[0] == mol.atom(0)
    
    with pytest.raises(KeyError):
        mol.atom("CA")

    s = mol.atoms("CA")

    assert len(s) == 351
    assert type(s) is not list
    assert s.what() == "SireMol::Selector<SireMol::Atom>"

    assert s[0].name().value() == "CA"
    assert mol[s[0].index()].index() == s[0].index()

    assert s == mol["atomname CA"]
    assert s == mol["CA"]
    assert s == mol[sr.atomid("CA")]

    s = mol.residues("resnum <= 10")

    assert len(s) == 10
    assert type(s) is not list
    assert s.what() == "SireMol::Selector<SireMol::Residue>"

    assert s[0].name().value() == "GLU"

    for i in range(0, 10):
        assert s[i].number().value() == i+1

    assert s == mol["resnum <= 10"]
    assert s == mol["resnum 1, 2, 3, 4, 5, 6, 7, 8, 9, 10"]

    assert mol["atoms in resnum 1"] == mol.residue(0).atoms()

    s = mol["residues with atomname HZ"]

    assert len(s) == 13

    for res in s:
        assert res.name().value() == "PHE"
        assert res["HZ"].name().value() == "HZ"
        assert res.atom("HZ").residue() == res

    from sire.mol import AtomName, AtomNum

    names = [AtomName('N'), AtomName('CA'), AtomName('C'), AtomName('O'), AtomName('CB')]

    assert mol[0:5].names() == names
    assert mol[4::-1].names() == names[4::-1]

    numbers = [AtomNum(7), AtomNum(9), AtomNum(11)]

    assert mol[6:12:2].numbers() == numbers
    assert mol[10:4:-2].numbers() == numbers[2::-1]

    with pytest.raises(KeyError):
        mol["X"]

    with pytest.raises(IndexError):
        mol[10000]

    s = mol["resnum >= 5 and resnum <= 10"]

    assert len(s) == 6

    for i in range(0, 6):
        assert s[i].number().value() == i+5

    s = mol[sr.resid(idx=0)]["atomname C, CA"]

    assert len(s) == 2

    assert s["CA"].name().value() == "CA"
    assert s["C"].name().value() == "C"

    count = 0

    for res in mol["residx < 3"]:
        for atom in res["atomname C, N"]:
            count += 1

    assert count == 6


def test_single_select_is_right(alanin_mols):
    import sire as sr
    mols = alanin_mols

    mol = mols[0]

    assert len(mol.segments()) == 1
    assert mol.segments().names()[0] == mol.segment(0).name()

    assert mol.segments("MAIN")[0].what() == sr.mol.Segment.typename()
    assert mol.segment("MAIN").what() == sr.mol.Segment.typename()
    assert mol["segname MAIN"].what() == sr.mol.Molecule.typename()


def test_sizes_correct(ala_mols):
    mols = ala_mols

    assert len(mols) == mols.count()
    assert len(mols) == mols.size()
    assert len(mols) == mols.num_molecules()

    for mol in mols:
        assert len(mol) == mol.count()
        assert len(mol) == mol.size()
        assert len(mol) == mol.num_atoms()

        for res in mol.residues():
            assert len(res) == res.count()
            assert len(res) == res.size()
            assert len(res) == res.num_atoms()


def test_search_terms(ala_mols):
    mols = ala_mols

    for atom in mols["mass >= 16 g_per_mol"]:
        assert atom.mass().value() >= 16.0

    for atom in mols["charge < 0"]:
        assert atom.charge().value() < 0

    for atom in mols["mass >= 16 and charge < 0"]:
        assert atom.charge().value() < 0
        assert atom.mass().value() >= 16

    check_mass = mols[0][0].mass().value()

    atoms = mols[0][f"mass >= {check_mass-0.001} and mass <= {check_mass+0.001}"]

    assert len(atoms) > 0

    for atom in atoms:
        assert atom.mass().value() == pytest.approx(check_mass)

    atoms = mols[0][f"mass =~ {check_mass}"]

    assert len(atoms) > 0

    for atom in atoms:
        assert atom.mass().value() == pytest.approx(check_mass)

    check_charge = mols[0][1].charge().value()

    atoms = mols[0][f"charge >= {check_charge-0.001} and charge <= {check_charge+0.001}"]

    assert len(atoms) > 0

    for atom in atoms:
        assert atom.charge().value() == pytest.approx(check_charge)

    atoms = mols[0][f"charge =~ {check_charge}"]

    assert len(atoms) > 0

    for atom in atoms:
        assert atom.charge().value() == pytest.approx(check_charge)

    import sire.search
    old_eps = sire.search.get_approx_epsilon()

    sire.search.set_approx_epsilon(1e-20)

    with pytest.raises(KeyError):
        mols[0][f"charge =~ {check_charge}"]

    sire.search.set_approx_epsilon(old_eps)

    assert sire.search.get_approx_epsilon() == old_eps


def test_in_searches(ala_mols):
    mols = ala_mols

    assert len(mols["atoms in *"]) == mols.num_atoms()
    assert len(mols["residues in *"]) == mols.num_residues()

    assert len(mols["atoms in resname ACE"]) == mols["resname ACE"].num_atoms()
    assert len(mols["atoms in (molecules with count(element O) > 1)"]) == mols[0].num_atoms()
    assert len(mols["atoms in (molecules with resname ALA)"]) == mols[0].num_atoms()
    assert len(mols["(atoms in molecules) with resname ALA"]) == mols["resname ALA"].num_atoms()

    assert mols["atoms in molidx 0"] == mols[0]["atoms"]
    assert mols["atoms in molidx -1"] == mols[-1]["atoms"]

    # check precedence - want "atoms in molecules with resname ALA" to 
    # be "atoms in (molecules with resname ALA)"
    assert mols["atoms in molecules with resname ALA"] == mols["atoms in (molecules with resname ALA)"]


def test_with_searches(ala_mols):
    mols = ala_mols

    import sire as sr

    for mol in mols["molecules with count(atoms) >= 3"]:
        assert(mol.num_atoms() >= 3)

    for res in mols["residues with count(atoms) >= 3"]:
        assert(res.num_atoms() >= 3)

    # single residue match
    res = mols["residues with count(atoms) > 6"]

    assert res.name().value() == "ALA"
    assert res.num_atoms() > 6

    for atom in mols["atoms with resname ALA"]:
        assert atom.residue().name().value() == "ALA"

    assert len(mols["residues with element C"]) == 3

    # find the one molecule with two oxygens
    assert len(mols["molecules with count(element O) > 1"]["element O"]) == 2

    # there are no residues that contain more than one oxygen
    with pytest.raises(KeyError):
        mols["residues with count(element O) > 1"]

    # all but one molecule contains at one oxygen
    assert len(mols["molecules with count(element O) == 1"]) == len(mols)-1

    assert mols["atoms with molidx 0"] == mols[0]["atoms"]
    assert mols["atoms with molidx -1"] == mols[-1]["atoms"]


def test_all_searches(ala_mols):
    mols = ala_mols

    assert len(mols["atoms"]) == 1912
    assert len(mols["residues"]) == 633

    assert len(mols["atoms in *"]) == 1912
    assert len(mols["residues in *"]) == 633

    assert len(mols["atoms in molidx 0"]) == 22
    assert len(mols["residues in molidx 0"]) == 3


def test_count_searches(ala_mols):
    mols = ala_mols

    for mol in mols["molecules with count(atoms) == 3"]:
        assert mol.num_atoms() == 3

    mol = mols["molecules with count(residues) == 3"]
    assert mol.num_residues() == 3

