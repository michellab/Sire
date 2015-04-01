
import Sire.Stream

from Sire.MM import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Qt import *
from Sire.Maths import *
from Sire.Units import *
from Sire.IO import *

from nose.tools import assert_almost_equal

import sys

realcompare = False

try:
    realcompare = int(sys.argv[1])
except:
    pass

if realcompare:
    amber = Amber();
    amber.set14Factors(0, 0)

    (molecules, space) = amber.readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    print(space)

    protein = molecules[ MolNum(1) ].molecule()

    print("Protein has %d atoms" % protein.nAtoms())
else:
    protein = Sire.Stream.load("../io/protein.s3")
    space = PeriodicBox( Vector(77.3667, 84.0572, 86.8795) )

    cnrg_vacuum = -16406.191630892663  
    ljnrg_vacuum = -3171.065732415715
    ljnrg_vacuum_geo = -3182.4573192579305

    cnrg_box = -16395.045120818653  
    ljnrg_box = -3142.701735153256
    ljnrg_box_geo = -3154.5987213176636

    cnrg_group = -81.68430463380852  
    ljnrg_group = -5.71092316497434
    ljnrg_group_geo = -5.576781489836948

    cnrg_group_box = -113.79643732412087
    ljnrg_group_box = -6.793124009768353
    ljnrg_group_box_geo = -6.63391491344806

    cnrg_space = cnrg_vacuum
    ljnrg_space = ljnrg_vacuum

space = PeriodicBox( Vector(58.5,58.5,58.5) )
group_space = PeriodicBox( Vector(45,45,45) )

coul_cutoff = 15 * angstrom
lj_cutoff = 15 * angstrom

intraff = IntraCLJFF("old")
intraff.setSpace( Cartesian() )
intraff.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,lj_cutoff) )
intraff.setShiftElectrostatics( True )

intraff.add(protein)

intraff01 = IntraGroupCLJFF("old01")
intraff01.setSpace( Cartesian() )
intraff01.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,lj_cutoff) )
intraff01.setShiftElectrostatics( True )

nres = 1

select0 = protein.selection()
select0 = select0.selectNone()
select1 = protein.selection()
select1 = select1.selectNone()

for i in range(0,nres):
    select0 = select0.select( ResIdx(i) )

for i in range(nres, protein.nResidues()):
    select1 = select1.select( ResIdx(i) )

group0 = PartialMolecule( protein, select0 )
group1 = PartialMolecule( protein, select1 )

intraff01.add( group0, MGIdx(0) )
intraff01.add( group1, MGIdx(1) )

intrafunc = CLJIntraShiftFunction(coul_cutoff, lj_cutoff)
intrafunc.setConnectivity(protein)
intrafunc.setSpace( Cartesian() )

cljatoms = CLJAtoms(protein, CLJAtoms.USE_ATOMIDX)
cljboxes = CLJBoxes(cljatoms, 7.5*angstrom)

cljatoms0 = CLJAtoms(group0, CLJAtoms.USE_ATOMIDX)
cljatoms1 = CLJAtoms(group1, CLJAtoms.USE_ATOMIDX)

cljboxes0 = CLJBoxes(cljatoms0, 7.5*angstrom)
cljboxes1 = CLJBoxes(cljatoms1, 7.5*angstrom)

def pvt_compare_group(verbose):

    t = QElapsedTimer()
    t.start()
    intraff01.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = intraff01.energy( intraff01.components().coulomb() ).value()
    old_ljnrg = intraff01.energy( intraff01.components().lj() ).value()

    if not realcompare:
        # here is the energy of the protein without the 1-4 terms
        # (the new forcefield doesn't calculate 1-4 terms)
        old_cnrg = cnrg_space
        old_ljnrg = ljnrg_space

    t.start()
    (new_cnrg, new_ljnrg) = intrafunc.calculate(cljatoms0, cljatoms1)
    new_ns = t.nsecsElapsed()

    cljcalculator = CLJCalculator()
    t.start()
    (box_cnrg, box_ljnrg) = cljcalculator.calculate(intrafunc, cljboxes0, cljboxes1)
    box_ns = t.nsecsElapsed()                

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))

    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )
    

def pvt_compare(verbose):

    t = QElapsedTimer()
    t.start()
    intraff.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = intraff.energy( intraff.components().coulomb() ).value()
    old_ljnrg = intraff.energy( intraff.components().lj() ).value()

    if not realcompare:
        # here is the energy of the protein without the 1-4 terms
        # (the new forcefield doesn't calculate 1-4 terms)
        old_cnrg = cnrg_space  
        old_ljnrg = ljnrg_space

    t.start()
    (new_cnrg, new_ljnrg) = intrafunc.calculate(cljatoms)
    new_ns = t.nsecsElapsed()

    cljcalculator = CLJCalculator()
    t.start()
    (box_cnrg, box_ljnrg) = cljcalculator.calculate(intrafunc, cljboxes)
    box_ns = t.nsecsElapsed()

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))
    
    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )


def test_compare_vacuum(verbose = False):
    intraff.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_vacuum
        globals()["ljnrg_space"] = ljnrg_vacuum

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare(verbose)

def test_compare_vacuum_geo(verbose = False):
    intraff.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )
    
    if not realcompare:
        globals()["cnrg_space"] = cnrg_vacuum
        globals()["ljnrg_space"] = ljnrg_vacuum_geo
    
    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )
    
    pvt_compare(verbose)

def test_compare_box(verbose = False):
    intraff.setSpace( space )
    intrafunc.setSpace( space )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_box
        globals()["ljnrg_space"] = ljnrg_box

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare(verbose)

def test_compare_box_geo(verbose = False):
    intraff.setSpace( space )
    intrafunc.setSpace( space )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_box
        globals()["ljnrg_space"] = ljnrg_box_geo

    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare(verbose)

def test_compare_group(verbose = False):
    intraff01.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_group
        globals()["ljnrg_space"] = ljnrg_group

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare_group(verbose)

def test_compare_group_geo(verbose = False):
    intraff01.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_group
        globals()["ljnrg_space"] = ljnrg_group_geo

    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare_group(verbose)

def test_compare_group_box(verbose = False):
    intraff01.setSpace( group_space )
    intrafunc.setSpace( group_space )

    if not realcompare:
        globals()["cnrg_space"] = cnrg_group_box
        globals()["ljnrg_space"] = ljnrg_group_box

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare_group(verbose)

def test_compare_group_box_geo(verbose = False):
    intraff01.setSpace( group_space )
    intrafunc.setSpace( group_space )
        
    if not realcompare:
        globals()["cnrg_space"] = cnrg_group_box
        globals()["ljnrg_space"] = ljnrg_group_box_geo
    
    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare_group(verbose)

if __name__ == "__main__":
    print("\nArithmetic combining rules")

    print("\nIntramolecular vacuum energy")
    test_compare_vacuum(True)

    print("\nIntramolecular periodic boundaries energy")
    test_compare_box(True)

    print("\nGroup intramolecular vacuum energy")
    test_compare_group(True)

    print("\nGroup intramolecular periodic boundaries energy")
    test_compare_group_box(True)

    print("\nGeometric combining rules")

    print("\nIntramolecular vacuum energy")
    test_compare_vacuum_geo(True)

    print("\nIntramolecular periodic boundaries energy")
    test_compare_box_geo(True)
        
    print("\nGroup intramolecular vacuum energy")
    test_compare_group_geo(True)
    
    print("\nGroup intramolecular periodic boundaries energy")
    test_compare_group_box_geo(True)

