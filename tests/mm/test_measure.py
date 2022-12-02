
import pytest


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


def test_periodic_measure(ala_mols):
    """Test that distances are calculated correctly taking
       into account periodic boundary conditions
    """
    import sire as sr
    import random

    mols = ala_mols

    space = mols.property("space")

    for i in range(0,100):
        idx0 = random.randint(0, len(mols)-1)
        idx1 = random.randint(0, len(mols)-1)

        c0 = mols[idx0].evaluate().center_of_mass()
        c1 = mols[idx1].evaluate().center_of_mass()

        # test that the periodic distance is calculated correctly
        dist0 = sr.measure(mols[idx0], mols[idx1])
        dist1 = space.calc_dist(c0, c1)

        assert dist0.value() == pytest.approx(dist1)

        # now test that the infinite cartesian distance is
        # calculated correctly
        dist0 = sr.measure(mols[idx0], mols[idx1], ignore_space=True)
        dist1 = sr.maths.Vector.distance(c0, c1)

        assert dist0.value() == pytest.approx(dist1.value())
