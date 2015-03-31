
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
exp_system = System()

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
        else:
            waters.add(water)

grid_system.add(swapwaters)
grid_system.add(waters)
exp_system.add(swapwaters)
exp_system.add(waters)

gridff = GridFF("gridff")
gridff.setBuffer(2 * angstrom)
gridff.setGridSpacing( 0.75 * angstrom )
gridff.setLJCutoff( 10 * angstrom )
gridff.setCoulombCutoff( 25 * angstrom )

gridff.add(swapwaters, MGIdx(0))
gridff.add(waters, MGIdx(1))
gridff.setSpace( Cartesian() )
gridff.setShiftElectrostatics(True)

cljff = InterGroupCLJFF("cljff")
cljff.setSwitchingFunction( HarmonicSwitchingFunction(25*angstrom, 25*angstrom, 10*angstrom, 10*angstrom) )
cljff.add(swapwaters, MGIdx(0))
cljff.add(waters, MGIdx(1))
cljff.setSpace( Cartesian() )
cljff.setShiftElectrostatics(True)

swap_swapff = InterCLJFF("swap-swap")
swap_swapff.setSpace(Cartesian())
swap_swapff.setSwitchingFunction( HarmonicSwitchingFunction(25*angstrom, 25*angstrom, 10*angstrom, 10*angstrom) )
swap_swapff.add(swapwaters)
grid_system.add(swap_swapff)
exp_system.add(swap_swapff)

grid_system.add(gridff)
exp_system.add(cljff)

grid_system.setComponent( grid_system.totalComponent(),  \
                          gridff.components().total() + swap_swapff.components().total() )
exp_system.setComponent( exp_system.totalComponent(), \
                              cljff.components().total() + swap_swapff.components().total() )

print((grid_system.energies()))
print((exp_system.energies()))

print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

diff = grid_system.energy() - exp_system.energy()
print(("The difference is %s\n" % diff))

rbmc = RigidBodyMC(swapwaters)
rbmc.setReflectionSphere(center_point, 7.5*angstrom)

moves = SameMoves(rbmc)

PDB().write(grid_system.molecules(), "test0000.pdb")

t = QTime()

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system = moves.move(grid_system, 1000, False)
    ms = t.elapsed()
    print(("Moves complete! Took %d ms" % ms))
    print(("GRID: ",grid_system.energies()))
    exp_system.update( grid_system.molecules() )
    print(("EXPT: ",exp_system.energies()))

    print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

    diff = grid_system.energy() - exp_system.energy()
    print(("The difference is %s\n" % diff))

    PDB().write(grid_system.molecules(), "test%0004d.pdb" % i)

# Save and restore the two systems from binary
import Sire.Stream
print("Saving the grid system...")
Sire.Stream.save( (grid_system, exp_system), "test/SireMM/testgrid.s3" )
print("Reloading the grid system...")
(grid_system, exp_system) = Sire.Stream.load("test/SireMM/testgrid.s3")


print((grid_system.energies()))
print((exp_system.energies()))
print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

diff = grid_system.energy() - exp_system.energy()
print(("The difference is %s\n" % diff))

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system = moves.move(grid_system, 1000, False)
    ms = t.elapsed()
    print(("Moves complete! Took %d ms" % ms))
    print(("GRID: ",grid_system.energies()))
    exp_system.update( grid_system.molecules() )
    print(("EXPT: ",exp_system.energies()))

    print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

    diff = grid_system.energy() - exp_system.energy()
    print(("The difference is %s\n" % diff))

    PDB().write(grid_system.molecules(), "test%0004d.pdb" % (i+10))

print("\nTriggering recalculation from scratch...")
grid_system.mustNowRecalculateFromScratch()
exp_system.mustNowRecalculateFromScratch()

print((grid_system.energies()))
print((exp_system.energies()))
print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

diff = grid_system.energy() - exp_system.energy()
print(("The difference is %s\n" % diff))

for i in range(1,11):
    print("Moving the system...")
    t.start()
    grid_system = moves.move(grid_system, 1000, False)
    ms = t.elapsed()
    print(("Moves complete! Took %d ms" % ms))
    print(("GRID: ",grid_system.energies()))
    exp_system.update( grid_system.molecules() )
    print(("EXPT: ",exp_system.energies()))

    print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))

    diff = grid_system.energy() - exp_system.energy()
    print(("The difference is %s\n" % diff))
    
    PDB().write(grid_system.molecules(), "test%0004d.pdb" % (i+20))

print("\nTriggering recalculation from scratch...")
grid_system.mustNowRecalculateFromScratch()
exp_system.mustNowRecalculateFromScratch()
    
print((grid_system.energies()))
print((exp_system.energies()))
print(("\nGrid energy equals: %s. Explicit energy equals: %s." % \
          (grid_system.energy(), exp_system.energy())))
    
diff = grid_system.energy() - exp_system.energy()
print(("The difference is %s\n" % diff))

print("Done :-)")
