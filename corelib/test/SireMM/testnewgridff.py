
from Sire.IO import *
from Sire.MM import *
from Sire.System import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.FF import *
from Sire.Move import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Qt import *

import os

amber = Amber()

(molecules, space) = amber.readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")

grid_system = System()
grid_system2 = System()

swapwaters = MoleculeGroup("swapwaters")
waters = MoleculeGroup("waters")

molnums = molecules.molNums();

for molnum in molnums:
    water = molecules[molnum].molecule()

    if water.residue().number() == ResNum(2025):
        center_water = water

swapwaters.add(center_water)

center_point = center_water.evaluate().center()

for molnum in molnums:
    if molnum != center_water.number():
        water = molecules[molnum].molecule()

        if Vector.distance(center_point, water.evaluate().center()) < 7.5:
            water = water.residue().edit().setProperty("PDB-residue-name", "SWP").commit()
            swapwaters.add(water)
#            if swapwaters.isEmpty():
#                swapwaters.add(water)
        else:
             waters.add(water)
#            if waters.isEmpty():
#                waters.add(water)

grid_system.add(swapwaters)
grid_system.add(waters)
grid_system2.add(swapwaters)
grid_system2.add(waters)

gridff = GridFF("gridff")
gridff.setBuffer(2 * angstrom)
gridff.setGridSpacing( 0.25 * angstrom )
gridff.setLJCutoff( 10 * angstrom )
gridff.setCoulombCutoff( 15 * angstrom )

gridff.add(swapwaters, MGIdx(0))
gridff.add(waters, MGIdx(1))
gridff.setSpace( Cartesian() )
gridff.setShiftElectrostatics(True)
#gridff.setUseReactionField(True)
#gridff.setUseAtomisticCutoff(True)

gridff2 = GridFF2("gridff2")
gridff2.setBuffer(2 * angstrom)
gridff2.setGridSpacing( 0.25 * angstrom )
gridff2.setLJCutoff( 10 * angstrom )
gridff2.setCoulombCutoff( 15 * angstrom )

gridff2.add(swapwaters, MGIdx(0))
gridff2.add(waters, MGIdx(1))
gridff2.setSpace( Cartesian() )
gridff2.setShiftElectrostatics(True)
#gridff2.setUseReactionField(True)
#gridff2.setUseAtomisticCutoff(True)

swap_swapff = InterCLJFF("swap-swap")
swap_swapff.setSpace(Cartesian())
swap_swapff.setSwitchingFunction( HarmonicSwitchingFunction(25*angstrom, 25*angstrom, 10*angstrom, 10*angstrom) )
swap_swapff.add(swapwaters)
grid_system.add(swap_swapff)
grid_system2.add(swap_swapff)

grid_system.add(gridff)
grid_system2.add(gridff2)

grid_system.setComponent( grid_system.totalComponent(),  \
                          gridff.components().total() + swap_swapff.components().total() )
grid_system2.setComponent( grid_system2.totalComponent(), \
                           gridff2.components().total() + swap_swapff.components().total() )

print(grid_system.energies())
print(grid_system2.energies())

print("\nOld Grid energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

diff = grid_system.energy() - grid_system2.energy()
print("The difference is %s\n" % diff)

rbmc = RigidBodyMC(swapwaters)
rbmc.setReflectionSphere(center_point, 7.5*angstrom)

moves = SameMoves(rbmc)

PDB().write(grid_system2.molecules(), "test0000.pdb")

t = QTime()

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system2 = moves.move(grid_system2, 1000, False)
    ms = t.elapsed()
    print("Moves complete! Took %d ms" % ms)
    print("NEW GRIDFF: ",grid_system2.energies())
    grid_system.update( grid_system2.molecules() )
    print("OLD GRIDFF: ",grid_system.energies())

    print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

    diff = grid_system.energy() - grid_system2.energy()
    print("The difference is %s\n" % diff)

    PDB().write(grid_system2.molecules(), "test%0004d.pdb" % i)

# Save and restore the two systems from binary
import Sire.Stream
print("Saving the systems...")
Sire.Stream.save( (grid_system, grid_system2), "test/SireMM/testgrid.s3" )
print("Reloading the grid system...")
(grid_system, grid_system2) = Sire.Stream.load("test/SireMM/testgrid.s3")


print(grid_system.energies())
print(grid_system2.energies())
print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

diff = grid_system.energy() - grid_system2.energy()
print("The difference is %s\n" % diff)

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system2 = moves.move(grid_system2, 1000, False)
    ms = t.elapsed()
    print("Moves complete! Took %d ms" % ms)
    print("NEW GRIDFF: ",grid_system2.energies())
    grid_system.update( grid_system2.molecules() )
    print("OLD GRIDFF: ",grid_system.energies())

    print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

    diff = grid_system.energy() - grid_system2.energy()
    print("The difference is %s\n" % diff)

    PDB().write(grid_system2.molecules(), "test%0004d.pdb" % (i+10))

print("\nTriggering recalculation from scratch...")
grid_system.mustNowRecalculateFromScratch()
grid_system2.mustNowRecalculateFromScratch()

print(grid_system.energies())
print(grid_system2.energies())
print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

diff = grid_system.energy() - grid_system2.energy()
print("The difference is %s\n" % diff)

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system2 = moves.move(grid_system2, 1000, False)
    ms = t.elapsed()
    print("Moves complete! Took %d ms" % ms)
    print("NEW GRIDFF: ",grid_system2.energies())
    grid_system.update( grid_system2.molecules() )
    print("OLD GRIDFF: ",grid_system.energies())

    print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))

    diff = grid_system.energy() - grid_system2.energy()
    print("The difference is %s\n" % diff)
    
    PDB().write(grid_system2.molecules(), "test%0004d.pdb" % (i+20))

print("\nTriggering recalculation from scratch...")
grid_system.mustNowRecalculateFromScratch()
grid_system2.mustNowRecalculateFromScratch()
    
print(grid_system.energies())
print(grid_system2.energies())
print("\nOld GridFF energy equals: %s. New GridFF energy equals: %s." % \
          (grid_system.energy(), grid_system2.energy()))
    
diff = grid_system.energy() - grid_system2.energy()
print("The difference is %s\n" % diff)

print("Done :-)")
