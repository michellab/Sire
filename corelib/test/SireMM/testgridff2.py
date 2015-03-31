
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

coul_cutoff = 20 * angstrom
lj_cutoff = 10 * angstrom

amber = Amber()

(molecules, space) = amber.readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")

system = System()

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

system.add(swapwaters)
system.add(waters)

gridff = GridFF("gridff")
gridff.setCombiningRules("arithmetic")
print("Combining rules are %s" % gridff.combiningRules())
gridff.setBuffer(2 * angstrom)
gridff.setGridSpacing( 0.5 * angstrom )
gridff.setLJCutoff(lj_cutoff)
gridff.setCoulombCutoff(coul_cutoff)
gridff.setShiftElectrostatics(True)
#gridff.setUseAtomisticCutoff(True)
#gridff.setUseReactionField(True)

cljgridff = CLJGrid()
cljgridff.setCLJFunction( CLJShiftFunction(coul_cutoff,lj_cutoff) )
cljgridff.setFixedAtoms( CLJAtoms(waters.molecules()) )
cljatoms = CLJAtoms(swapwaters.molecules())
cljgridff.setGridDimensions( cljatoms, 0.5 * angstrom, 2 * angstrom )
print("Grid box equals %s" % cljgridff.grid())
cljboxes = CLJBoxes(cljatoms)
(cnrg, ljnrg) = cljgridff.calculate(cljboxes)
print("CLJGridFF:  %s  %s  %s" % (cnrg+ljnrg, cnrg, ljnrg))
cljgridff.setUseGrid(False)
(cnrg, ljnrg) = cljgridff.calculate(cljboxes)
print("CLJGridFF:  %s  %s  %s" % (cnrg+ljnrg, cnrg, ljnrg))


gridff.add(swapwaters, MGIdx(0))
gridff.add(waters, MGIdx(1))
gridff.setSpace( Cartesian() )

gridff2 = GridFF2("gridff2")
gridff2.setCombiningRules("arithmetic")
gridff2.setBuffer(2*angstrom)
gridff2.setGridSpacing( 0.5 * angstrom )
gridff2.setLJCutoff(lj_cutoff)
gridff2.setCoulombCutoff(coul_cutoff)
gridff2.setShiftElectrostatics(True)
#gridff2.setUseAtomisticCutoff(True)
#gridff2.setUseReactionField(True)

gridff2.add( swapwaters, MGIdx(0) )
gridff2.addFixedAtoms(waters.molecules())
gridff2.setSpace( Cartesian() )

testff = TestFF()
testff.add( swapwaters.molecules() )
testff.addFixedAtoms(waters.molecules())
testff.setCutoff(coul_cutoff, lj_cutoff)

cljff = InterGroupCLJFF("cljff")
cljff.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff) )
cljff.add(swapwaters, MGIdx(0))
cljff.add(waters, MGIdx(1))
cljff.setShiftElectrostatics(True)
#cljff.setUseAtomisticCutoff(True)
#cljff.setUseReactionField(True)
cljff.setSpace( Cartesian() )

cljff2 = InterCLJFF("cljff2")
cljff2.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff) )
cljff2.add(waters)
cljff2.setShiftElectrostatics(True)
cljff2.setSpace( Cartesian() )

print(gridff.energies())
print(gridff2.energies())

print("\nEnergies")
print(gridff.energies())
print(gridff2.energies())

t = QTime()
t.start()
nrgs = cljff.energies()
ms = t.elapsed()
print(cljff.energies())
print("Took %d ms" % ms)

testff.calculateEnergy()

t.start()
nrgs = cljff2.energies()
ms = t.elapsed()

print("\nExact compare")
print(cljff2.energies())
print("Took %d ms" % ms)
