
from Sire.MM import *
from Sire.FF import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Base import *
from Sire.Units import *
from Sire.Qt import *

from nose.tools import assert_almost_equal

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

grid_spacing = 0.5 * angstrom
grid_buffer = 2.0 * angstrom

nmoves = 1000

#space = Cartesian()

reflect_sphere_center = Vector(10)
reflect_sphere_radius = 10 * angstrom

long_coul_cutoff = 25 * angstrom
coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)
grid_switchfunc = HarmonicSwitchingFunction(long_coul_cutoff,long_coul_cutoff,
                                            lj_cutoff,lj_cutoff)

cluster = MoleculeGroup("cluster")
waters = MoleculeGroup("waters")

for i in range(0, mols.nMolecules()):
    mol = mols[ MolIdx(i) ].molecule()

    if Vector.distance(mol.evaluate().center(),reflect_sphere_center) < reflect_sphere_radius.value() :
        cluster.add(mol)
    else:
        waters.add(mol)

oldff = InterCLJFF("oldff")
oldff.setSwitchingFunction(switchfunc)
oldff.setSpace(space)
oldff.setShiftElectrostatics(True)

newff = InterFF("newff")
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)

old_clusterff = InterCLJFF("old_clusterff")
old_clusterff.setSwitchingFunction(grid_switchfunc)
old_clusterff.setSpace(Cartesian())
old_clusterff.setShiftElectrostatics(True)
old_fixedff = GridFF("old_fixedff")
old_fixedff.setSwitchingFunction(grid_switchfunc)
old_fixedff.setCoulombCutoff(long_coul_cutoff)
old_fixedff.setLJCutoff(lj_cutoff)
old_fixedff.setSpace(Cartesian())
old_fixedff.setShiftElectrostatics(True)
old_fixedff.setGridSpacing(grid_spacing)
old_fixedff.setBuffer(grid_buffer)

new_clusterff = InterFF("new_clusterff")
new_clusterff.setProperty("cljFunction", CLJShiftFunction())
new_clusterff.setProperty("switchingFunction", grid_switchfunc )
new_clusterff.setProperty("space", Cartesian())
new_clusterff.setProperty("gridSpacing", LengthProperty(grid_spacing))
new_clusterff.setProperty("gridBuffer", LengthProperty(grid_buffer))

old_clusterff.add(cluster)
old_fixedff.add(cluster, MGIdx(0))
old_fixedff.addFixedAtoms(waters)

new_clusterff.add(cluster)
new_clusterff.setFixedAtoms(waters.molecules())
new_clusterff.setUseGrid(True)

t = QElapsedTimer()
t.start()
cljmols = CLJBoxes( CLJAtoms(mols.molecules()) )
cljns = t.nsecsElapsed()

t.start()
oldff.add(mols)
oldns = t.nsecsElapsed()

t.start()
newff.add(mols)
newns = t.nsecsElapsed()

print("Setup times: %s ms vs. %s ms (%s ms)" % (0.000001*oldns,0.000001*newns,
                                                0.000001*cljns))

def test_energy(verbose = False):
    t = QElapsedTimer()

    cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    cljfunc.setSpace(space)
    cljcalc = CLJCalculator()
    t.start()
    (cnrg,ljnrg) = cljcalc.calculate(cljfunc,cljmols)
    ns = t.nsecsElapsed()

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()
    oldcnrg = oldff.energy(oldff.components().coulomb()).value()
    oldljnrg = oldff.energy(oldff.components().lj()).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy(newff.components().coulomb()).value()
    newljnrg = newff.energy(newff.components().lj()).value()

    if verbose:
        print("\nTotal energy")
        print("CLJFUNC:  %s  %s  %s  : %s ms" % (cnrg+ljnrg,cnrg,ljnrg,
                                                 0.000001*ns))
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))

    assert_almost_equal( cnrg, newcnrg, 6 )
    assert_almost_equal( ljnrg, newljnrg, 6 )
    assert_almost_equal( oldcnrg, newcnrg, 2 )
    assert_almost_equal( oldljnrg, newljnrg, 2 )

    water = mols[ MolIdx(0) ].molecule()
    water = water.move().translate( Vector(1,0,0) ).commit()

    t.start()
    oldff.update(water)
    old_change = t.nsecsElapsed()

    t.start()
    newff.update(water)
    new_change = t.nsecsElapsed()

    if verbose:
        print("ff.update(): %s ms vs. %s ms" % (old_change*0.000001,
                                                new_change*0.000001))

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()
    oldcnrg = oldff.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldff.energy( oldff.components().lj() ).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy( newff.components().coulomb() ).value()
    newljnrg = newff.energy( newff.components().lj() ).value()

    if verbose:
        print("\nChanged energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))

    assert_almost_equal( oldcnrg, newcnrg, 1 )
    assert_almost_equal( oldljnrg, newljnrg, 1 )

    oldff.mustNowRecalculateFromScratch()
    newff.mustNowRecalculateFromScratch()

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()

    r_oldcnrg = oldff.energy( oldff.components().coulomb() ).value()
    r_oldljnrg = oldff.energy( oldff.components().lj() ).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()

    r_newcnrg = newff.energy( newff.components().coulomb() ).value()
    r_newljnrg = newff.energy( newff.components().lj() ).value()

    if verbose:
        print("\nRecalculated energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (r_oldcnrg+r_oldljnrg,r_oldcnrg,r_oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))

    assert_almost_equal( oldcnrg, r_oldcnrg, 6 )
    assert_almost_equal( oldljnrg, r_oldljnrg, 6 )
    assert_almost_equal( newcnrg, r_newcnrg, 6 )
    assert_almost_equal( newljnrg, r_newljnrg, 6 )

def test_sim(verbose = False):

    oldsys = System()
    newsys = System()

    #oldsys.add(mols)
    #newsys.add(mols)

    oldsys.add(oldff)
    newsys.add(newff)

    t = QElapsedTimer()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()

    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()

    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()

    if verbose:
        print("\nStarting energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))

    moves = RigidBodyMC(mols)
    moves.setGenerator( RanGenerator( 42 ) )

    t.start()
    moves.move(oldsys, nmoves, False)
    move_oldns = t.nsecsElapsed()

    old_naccepted = moves.nAccepted()
    old_nrejected = moves.nRejected()

    moves.setGenerator( RanGenerator( 42 ) )
    moves.clearStatistics()

    t.start()
    moves.move(newsys, nmoves, False)
    move_newns = t.nsecsElapsed()

    new_naccepted = moves.nAccepted()
    new_nrejected = moves.nRejected()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()
    
    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()

    if verbose:
        print("\nMoves: %s ms vs. %s ms" % (0.000001*move_oldns, 0.000001*move_newns))
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("nAccepted() = %s, nRejected() = %s" % (old_naccepted, old_nrejected))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("nAccepted() = %s, nRejected() = %s" % (new_naccepted, new_nrejected))
    
    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    
    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    r_oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    r_oldljnrg = oldsys.energy( oldff.components().lj() ).value()

    r_newcnrg = newsys.energy( newff.components().coulomb() ).value()
    r_newljnrg = newsys.energy( newff.components().lj() ).value()

    if verbose:
        print("\nRecalculated energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (r_oldcnrg+r_oldljnrg,r_oldcnrg,r_oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))

def test_fixed_sim(verbose = False):
    oldsys = System()
    newsys = System()

    oldsys.add(cluster)
    oldsys.add(old_fixedff)

    newsys.add(cluster)
    new_fixedff = new_clusterff.clone()
    new_fixedff.setProperty("fixedOnly", BooleanProperty(True))
    newsys.add(new_fixedff)

    t = QElapsedTimer()
    t.start()
    old_total = oldsys.energy().value()
    oldns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    newns = t.nsecsElapsed()

    ff = newsys[FFName("new_clusterff")]
    print(ff.grid())

    old_cnrg = oldsys.energy( old_fixedff.components().coulomb() ).value()
    old_ljnrg = oldsys.energy( old_fixedff.components().lj() ).value() 

    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("OLD:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                             0.000001*oldns))
        print("NEW:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                             0.000001*newns))

    moves = RigidBodyMC(cluster)                    
    moves.setReflectionSphere( reflect_sphere_center, reflect_sphere_radius )
    moves.setGenerator( RanGenerator( 42 ) )
    
    t.start()
    moves.move(oldsys, 1000, False)
    move_oldns = t.nsecsElapsed()

    moves.setGenerator( RanGenerator( 42 ) )

    t.start()
    moves.move(newsys, 1000, False)
    move_newns = t.nsecsElapsed()

    t.start()
    old_total = oldsys.energy().value()
    old_ns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    new_ns = t.nsecsElapsed()

    old_cnrg = oldsys.energy( old_fixedff.components().coulomb() ).value()    
    old_ljnrg = oldsys.energy( old_fixedff.components().lj() ).value()    

    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("\nMoves: %s ms vs. %s ms" % (0.000001*move_oldns, 0.000001*move_newns))
        print("OLD SYS:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                                 0.000001*old_ns))
        print("NEW SYS:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                                 0.000001*new_ns))

    newsys.mustNowRecalculateFromScratch()
    oldsys.mustNowRecalculateFromScratch()

    t.start()
    old_total = oldsys.energy().value()
    old_ns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    new_ns = t.nsecsElapsed()

    old_cnrg = oldsys.energy( old_fixedff.components().coulomb() ).value()    
    old_ljnrg = oldsys.energy( old_fixedff.components().lj() ).value()    
    
    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("\nRecalculate energy")
        print("OLD SYS:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                                 0.000001*old_ns))
        print("NEW SYS:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                                 0.000001*new_ns))


def test_grid_sim(verbose = False):
    oldsys = System()
    newsys = System()

    oldsys.add(cluster)
    oldsys.add(old_clusterff)
    oldsys.add(old_fixedff)

    newsys.add(cluster)
    newsys.add(new_clusterff)

    t = QElapsedTimer()
    t.start()
    old_total = oldsys.energy().value()
    oldns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    newns = t.nsecsElapsed()

    ff = newsys[FFName("new_clusterff")]
    print(ff.grid())

    old_cnrg = oldsys.energy( old_clusterff.components().coulomb() ).value() + \
               oldsys.energy( old_fixedff.components().coulomb() ).value()
    old_ljnrg = oldsys.energy( old_clusterff.components().lj() ).value() + \
                oldsys.energy( old_fixedff.components().lj() ).value()

    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("OLD:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                             0.000001*oldns))
        print("NEW:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                             0.000001*newns))

    moves = RigidBodyMC(cluster)                    
    moves.setReflectionSphere( reflect_sphere_center, reflect_sphere_radius )
    moves.setGenerator( RanGenerator( 42 ) )
    
    t.start()
    moves.move(oldsys, 1000, False)
    move_oldns = t.nsecsElapsed()

    moves.setGenerator( RanGenerator( 42 ) )

    t.start()
    moves.move(newsys, 1000, False)
    move_newns = t.nsecsElapsed()

    t.start()
    old_total = oldsys.energy().value()
    old_ns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    new_ns = t.nsecsElapsed()

    old_cnrg = oldsys.energy( old_clusterff.components().coulomb() ).value() + \
               oldsys.energy( old_fixedff.components().coulomb() ).value()
    old_ljnrg = oldsys.energy( old_clusterff.components().lj() ).value() + \
                oldsys.energy( old_fixedff.components().lj() ).value()

    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("\nMoves: %s ms vs. %s ms" % (0.000001*move_oldns, 0.000001*move_newns))
        print("OLD SYS:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                                 0.000001*old_ns))
        print("NEW SYS:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                                 0.000001*new_ns))

    newsys.mustNowRecalculateFromScratch()
    oldsys.mustNowRecalculateFromScratch()

    t.start()
    old_total = oldsys.energy().value()
    old_ns = t.nsecsElapsed()

    t.start()
    new_total = newsys.energy().value()
    new_ns = t.nsecsElapsed()

    old_cnrg = oldsys.energy( old_clusterff.components().coulomb() ).value() + \
               oldsys.energy( old_fixedff.components().coulomb() ).value()
    old_ljnrg = oldsys.energy( old_clusterff.components().lj() ).value() + \
                oldsys.energy( old_fixedff.components().lj() ).value()
    
    new_cnrg = newsys.energy( new_clusterff.components().coulomb() ).value()
    new_ljnrg = newsys.energy( new_clusterff.components().lj() ).value()

    if verbose:
        print("\nRecalculate energy")
        print("OLD SYS:  %s  %s  %s  %s  : %s ms" % (old_total,old_cnrg+old_ljnrg,old_cnrg,old_ljnrg,
                                                 0.000001*old_ns))
        print("NEW SYS:  %s  %s  %s  %s  : %s ms" % (new_total,new_cnrg+new_ljnrg,new_cnrg,new_ljnrg,
                                                 0.000001*new_ns))


if __name__ == "__main__":
    test_energy(True)
    test_sim(True)
    test_grid_sim(True)
    test_fixed_sim(True)
