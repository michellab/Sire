
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Qt import *

from nose.tools import assert_almost_equal

coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

compare_energies = False

if compare_energies:
    (waters, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top", "peratom")
else:
    (waters, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

group0 = MoleculeGroup("group0")
group1 = MoleculeGroup("group1")

nwaters = 1

for i in range(0,nwaters):
    group0.add( waters[MolIdx(i)] )

for i in range(nwaters, waters.nMolecules()):
    group1.add( waters[MolIdx(i)] )

cljatoms = CLJAtoms(waters.molecules())
cljatoms0 = CLJAtoms(group0.molecules())
cljatoms1 = CLJAtoms(group1.molecules())

boxes = CLJBoxes(cljatoms)
boxes0 = CLJBoxes(cljatoms0)
boxes1 = CLJBoxes(cljatoms1)

cljff = InterCLJFF("all")
cljff.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff) )
cljff.add(waters)
cljff.setShiftElectrostatics(True)

cljff01 = InterGroupCLJFF("pair")
cljff01.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff) )
cljff01.add(group0, MGIdx(0))
cljff01.add(group1, MGIdx(1))
cljff01.setShiftElectrostatics(True)

cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)

def pvt_compare_group(verbose):

    cljcalc = CLJCalculator()

    t = QElapsedTimer()
    t.start()
    cljff01.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = cljff01.energy(cljff01.components().coulomb()).value()
    old_ljnrg = cljff01.energy(cljff01.components().lj()).value()

    t.start()
    (new_cnrg, new_ljnrg) = cljfunc.calculate(cljatoms0, cljatoms1)
    new_ns = t.nsecsElapsed()

    t.start()
    (box_cnrg, box_ljnrg) = cljfunc.calculate(boxes0, boxes1)
    box_ns = t.nsecsElapsed()

    t.start()
    (cal_cnrg, cal_ljnrg) = cljcalc.calculate(cljfunc, boxes0, boxes1)
    cal_ns = t.nsecsElapsed()

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))
        print("CAL: %s  %s  %s  %s ms" % (cal_cnrg+cal_ljnrg, cal_cnrg, cal_ljnrg, 0.000001*cal_ns))

    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( cal_cnrg, old_cnrg, 2 )
    assert_almost_equal( cal_ljnrg, cal_ljnrg, 2 )


def pvt_compare(verbose):

    cljcalc = CLJCalculator()

    t = QElapsedTimer()
    t.start()
    cljff.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = cljff.energy(cljff.components().coulomb()).value()
    old_ljnrg = cljff.energy(cljff.components().lj()).value()

    t.start()
    (new_cnrg, new_ljnrg) = cljfunc.calculate(cljatoms)
    new_ns = t.nsecsElapsed()

    t.start()
    (box_cnrg, box_ljnrg) = cljfunc.calculate(boxes)
    box_ns = t.nsecsElapsed()

    t.start()
    (cal_cnrg, cal_ljnrg) = cljcalc.calculate(cljfunc, boxes)
    cal_ns = t.nsecsElapsed()

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))
        print("CAL: %s  %s  %s  %s ms" % (cal_cnrg+cal_ljnrg, cal_cnrg, cal_ljnrg, 0.000001*cal_ns))

    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( cal_cnrg, old_cnrg, 2 )
    assert_almost_equal( cal_ljnrg, old_ljnrg, 2 )


def test_compare_vacuum(verbose = True):
    cljff.setSpace(Cartesian())
    cljff01.setSpace(Cartesian())
    cljfunc.setSpace(Cartesian())
    cljff.setCombiningRules("arithmetic")
    cljff01.setCombiningRules("arithmetic")
    cljfunc.setCombiningRules( CLJFunction.ARITHMETIC )
    pvt_compare(verbose)

def test_compare_vacuum_geo(verbose = True):
    cljff.setSpace(Cartesian())
    cljff01.setSpace(Cartesian())
    cljfunc.setSpace(Cartesian())
    cljff.setCombiningRules("geometric")
    cljff01.setCombiningRules("geometric")
    cljfunc.setCombiningRules( CLJFunction.GEOMETRIC )
    pvt_compare(verbose)

def test_compare_box(verbose = True):
    cljff.setSpace(space)
    cljff01.setSpace(space)
    cljfunc.setSpace(space)
    cljff.setCombiningRules("arithmetic")
    cljff01.setCombiningRules("arithmetic")
    cljfunc.setCombiningRules( CLJFunction.ARITHMETIC )
    pvt_compare(verbose)

def test_compare_box_geo(verbose = True):
    cljff.setSpace(space)
    cljff01.setSpace(space)
    cljfunc.setSpace(space)
    cljff.setCombiningRules("geometric")
    cljff01.setCombiningRules("geometric")
    cljfunc.setCombiningRules( CLJFunction.GEOMETRIC )
    pvt_compare(verbose)
    
def test_compare_vacuum_group(verbose = True):
    cljff.setSpace(Cartesian())
    cljff01.setSpace(Cartesian())
    cljfunc.setSpace(Cartesian())
    cljff.setCombiningRules("arithmetic")
    cljff01.setCombiningRules("arithmetic")
    cljfunc.setCombiningRules( CLJFunction.ARITHMETIC )
    pvt_compare_group(verbose)

def test_compare_vacuum_group_geo(verbose = True):
    cljff.setSpace(Cartesian())
    cljff01.setSpace(Cartesian())
    cljfunc.setSpace(Cartesian())
    cljff.setCombiningRules("geometric")
    cljff01.setCombiningRules("geometric")
    cljfunc.setCombiningRules( CLJFunction.GEOMETRIC )
    pvt_compare_group(verbose)

def test_compare_box_group(verbose = True):
    cljff.setSpace(space)
    cljff01.setSpace(space)
    cljfunc.setSpace(space)
    cljff.setCombiningRules("arithmetic")
    cljff01.setCombiningRules("arithmetic")
    cljfunc.setCombiningRules( CLJFunction.ARITHMETIC )
    pvt_compare_group(verbose)
    
def test_compare_box_group_geo(verbose = True):
    cljff.setSpace(space)
    cljff01.setSpace(space)
    cljfunc.setSpace(space)
    cljff.setCombiningRules("geometric")
    cljff01.setCombiningRules("geometric")
    cljfunc.setCombiningRules( CLJFunction.GEOMETRIC )
    pvt_compare_group(verbose)

if __name__ == "__main__":
    print("nCutGroups() = %s" % waters[MolIdx(0)].molecule().nCutGroups())
    print("space = %s" % space)

    print("\nArithmetic combining rules")
    print("\nVacuum boundary conditions")
    test_compare_vacuum(True)

    print("\nPeriodic box boundary conditions")
    test_compare_box(True)

    print("\nGroup Vacuum boundary conditions")
    test_compare_vacuum_group(True)

    print("\nGroup Periodic box boundary conditions")
    test_compare_box_group(True)

    print("\nGeometric combining rules")

    print("\nVacuum boundary conditions")
    test_compare_vacuum_geo(True)

    print("\nPeriodic box boundary conditions")
    test_compare_box_geo(True)

    print("\nGroup Vacuum boundary conditions")
    test_compare_vacuum_group_geo(True)
    
    print("\nGroup Periodic box boundary conditions")
    test_compare_box_group_geo(True)

