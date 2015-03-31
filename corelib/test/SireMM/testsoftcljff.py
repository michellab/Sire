
from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.Base import *
from Sire.Units import *

print("Reading in the coordinates of the waterbox...")
amber = Amber()
(waters, space) = amber.readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")
print("...done!")

cljff = InterCLJFF("hard")
softcljff = InterSoftCLJFF("soft")

cljff.add(waters)
softcljff.add(waters)

def buildSystem(forcefields):
    system = System()
    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("switchingFunction", HarmonicSwitchingFunction(8*angstrom, 7.5*angstrom))

    system.setProperty("shiftDelta", VariantProperty(2.0))
    system.setProperty("coulombPower", VariantProperty(0))

    return system

def printEnergies(alpha, system, components):

    energies = system.energies()

    print("%8.3f " % alpha, end=' ')

    for component in components:
        print("%14.10f " % energies[component], end=' ')

    print("\n", end=' ')

def printHeader(components):
    print("Alpha    ", end=' ')
    for component in components:
        print("%14s " % str(component), end=' ')

    print("\n", end=' ')

components = [ cljff.components().total(), softcljff.components().total(),
               cljff.components().coulomb(), softcljff.components().coulomb(),
               cljff.components().lj(), softcljff.components().lj() ]

print("\nTesting the group-based cutoff code...")

system = buildSystem( [cljff, softcljff] )

printHeader(components)

for i in range(0,11):
    alpha = 0.1 * i
    system.setProperty("alpha", VariantProperty(alpha))

    printEnergies(alpha, system, components)

print("\nTesting the force-shifted cutoff...")
cljff.setShiftElectrostatics(True)
softcljff.setShiftElectrostatics(True)

system = buildSystem( [cljff, softcljff] )

for i in range(0,11):
    alpha = 0.1 * i
    system.setProperty("alpha", VariantProperty(alpha))

    printEnergies(alpha, system, components)

print("\nTesting the reaction field cutoff...")
cljff.setUseReactionField(True)
softcljff.setUseReactionField(True)
cljff.setReactionFieldDielectric(78.3)
softcljff.setReactionFieldDielectric(78.3)

system = buildSystem( [cljff, softcljff] )

for i in range(0,11):
    alpha = 0.1 * i
    system.setProperty("alpha", VariantProperty(alpha))
    
    printEnergies(alpha, system, components)

