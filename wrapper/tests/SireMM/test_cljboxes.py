
from Sire.MM import *
from Sire.IO import *
from Sire.Mol import *
from Sire.Units import *
from Sire.Qt import *

from nose.tools import assert_almost_equal

(waters,space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

coul_cutoff = 25 * angstrom
lj_cutoff = 10 * angstrom

cljwaters = []

for i in range(0,waters.nMolecules()):
    cljwaters.append( CLJAtoms( waters[MolIdx(i)] ) )


def test_boxing(verbose = False):
    nloops = 10
    n0 = 100
    n1 = 101

    water = cljwaters[n0]

    boxes = CLJBoxes()

    if verbose:
        print("\nAdding this water to a box multiple times...")
        print(water)

    t = QElapsedTimer()
    add_ns = 0

    for i in range(0,nloops):
        t.start()
        idxs = boxes.add(water)
        add_ns += t.nsecsElapsed()

        assert( len(idxs) == water.count() )

        if verbose:
            print("nAtoms() == %s, %s x %s" % (boxes.nAtoms(), i+1,
                                               water.nAtoms()))
            print("idxs == %s" % idxs[0:3])

        assert( boxes.nAtoms() == (i+1)*water.nAtoms() )

        for j in range(0, water.nAtoms()):
            atom = water[j]
            idx = idxs[j]
            boxatom = boxes.at(idx)

            if verbose:
                print("%s: %s == %s" % (idx, atom, boxatom))

            assert( 1 + idx.index() > i * water.nAtoms() )
            assert( 1 + idx.index() <= (i+1) * water.nAtoms() )

            assert( atom.coordinates() == boxatom.coordinates() )
            assert( atom.charge() == boxatom.charge() )
            assert( atom == boxatom )

        for j in range(water.nAtoms(), water.count()):
            assert( water[j].isDummy() )
            assert( idxs[j].isNull() )

    if verbose:
        print("\nAverage adding time: %s ms" % (0.000001 * add_ns / nloops))

    cljwater0 = CLJBoxes( cljwaters[n0] )
    cljwater1 = CLJBoxes( cljwaters[n1] )

    cljcalc = CLJCalculator()
    cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)

    t.start()
    (cnrg,ljnrg) = cljcalc.calculate(cljfunc, cljwater0, cljwater1)
    ns = t.nsecsElapsed()

    t.start()
    (ncnrg,nljnrg) = cljcalc.calculate(cljfunc, boxes, cljwater1 )
    nns = t.nsecsElapsed()

    if verbose:
        print("\nONE NRG:  %s  %s  %s  : %s ms" % (cnrg+ljnrg,cnrg,ljnrg,0.000001*ns))
        print("\nALL NRG:  %s  %s  %s  : %s ms" % (ncnrg+nljnrg,ncnrg,nljnrg,0.000001*nns))

    assert_almost_equal( ncnrg, 10 * cnrg, 4 )
    assert_almost_equal( nljnrg, 10 * ljnrg, 4 )

if __name__ == "__main__":
    test_boxing(True)

