
from Sire.MM import *
from Sire.Mol import *
from Sire.IO import *
from Sire.Qt import *
from Sire.Units import *

from nose.tools import assert_almost_equal

coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

(waters, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

water = waters[ MolIdx(0) ].molecule()

waters.remove(water)

cljwaters = CLJBoxes( CLJAtoms(waters) )

def test_cljnegate(verbose = False):
    pos_atoms = CLJAtoms(water)

    t = QElapsedTimer()
    t.start()
    neg_atoms = pos_atoms.negate()
    ns = t.nsecsElapsed()

    if verbose:
        print("\nNegation took %s ms" % (0.000001*ns))
    
    assert( pos_atoms.count() == neg_atoms.count() )
    
    for i in range(len(pos_atoms)):
        assert( pos_atoms[i].charge() == -(neg_atoms[i].charge()) )   

def test_negate(verbose = False):
    cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    cljfunc.setSpace(space)

    pos_atoms = CLJAtoms(water)
    neg_atoms = pos_atoms.negate()

    all_atoms = (pos_atoms + neg_atoms).squeeze()

    pos_boxes = CLJBoxes(pos_atoms)
    neg_boxes = CLJBoxes(neg_atoms)
    all_boxes = CLJBoxes(all_atoms)

    t = QElapsedTimer()

    cljcalc = CLJCalculator()

    t.start()
    (tot_cnrg,tot_ljnrg) = cljcalc.calculate(cljfunc, cljwaters)
    tot_ns = t.nsecsElapsed()

    t.start()
    (pos_cnrg,pos_ljnrg) = cljcalc.calculate(cljfunc, pos_boxes, cljwaters)
    pos_ns = t.nsecsElapsed()

    t.start()
    (neg_cnrg,neg_ljnrg) = cljcalc.calculate(cljfunc, neg_boxes, cljwaters)
    neg_ns = t.nsecsElapsed()

    t.start()
    (all_cnrg,all_ljnrg) = cljcalc.calculate(cljfunc, all_boxes, cljwaters)
    all_ns = t.nsecsElapsed()

    if verbose:
        print("\nTOT:  %s  %s  %s  : %s ms" % (tot_cnrg+tot_ljnrg,
                                               tot_cnrg, tot_ljnrg,
                                               (0.000001*tot_ns)))
        print("POS:  %s  %s  %s  : %s ms" % (pos_cnrg+pos_ljnrg,
                                             pos_cnrg, pos_ljnrg,
                                             (0.000001*pos_ns)))
        print("NEG:  %s  %s  %s  : %s ms" % (neg_cnrg+neg_ljnrg,
                                             neg_cnrg, neg_ljnrg,
                                             (0.000001*neg_ns)))
        print("ALL:  %s  %s  %s  : %s ms" % (all_cnrg+all_ljnrg,
                                             all_cnrg, all_ljnrg,
                                             (0.000001*all_ns)))

    assert_almost_equal( pos_cnrg, -neg_cnrg, 6 )
    assert_almost_equal( pos_ljnrg, -neg_ljnrg, 6 )
    assert_almost_equal( all_cnrg, 0.0, 6 )
    assert_almost_equal( all_ljnrg, 0.0, 6 )

if __name__ == "__main__":
    test_cljnegate(True)
    test_negate(True)
