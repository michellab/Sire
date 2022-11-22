

def test_measure(ala_mols):
    import sire as sr

    mols = ala_mols

    bond = mols[0].bonds()[0]

    assert bond.measure() == sr.measure(bond[0], bond[1])

    ang = mols[0].angles()[0]

    assert ang.measure() == sr.measure(ang[0], ang[1], ang[2])

    dih = mols[0].dihedrals()[0]

    assert dih.measure() == sr.measure(dih[0], dih[1], dih[2], dih[3])

    assert dih.measure() == sr.measure(dih[0], dih[1], dih[2], dih[3],
                                       improper_angle=False)

    assert dih.measure() != sr.measure(dih[0], dih[1], dih[2], dih[3],
                                       improper_angle=True)

    imp = mols[0].impropers()[0]

    assert imp.measure() == sr.measure(imp[0], imp[1], imp[2], imp[3],
                                       improper_angle=False)

    assert imp.phi() == sr.measure(imp[0], imp[1], imp[2], imp[3],
                                   improper_angle=False)

    assert imp.theta() == sr.measure(imp[0], imp[1], imp[2], imp[3],
                                     improper_angle=True)
