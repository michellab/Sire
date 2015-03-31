
from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.Units import *

amber = Amber()

print("Loading the molecules...")
(molecules, space) = amber.readCrdTop("test/io/SYSTEM.crd", "test/io/SYSTEM.top")
print("...all loaded :-)")

protein = molecules[MolNum(2)].molecule()
#protein = molecules[MolNum(1)].molecule()

group_clj = IntraCLJFF("group_clj")

shift_clj = IntraCLJFF("shift_clj")
shift_clj.setShiftElectrostatics(True)

field_clj = IntraCLJFF("field_clj")
field_clj.setUseReactionField(True)
field_clj.setReactionFieldDielectric(78.3)

forcefields = [ group_clj, shift_clj, field_clj ]

group_coul = group_clj.components().coulomb()
shift_coul = shift_clj.components().coulomb()
field_coul = field_clj.components().coulomb()

system = System()

for forcefield in forcefields:
    forcefield.add(protein)
    system.add(forcefield)

def printEnergies(nrgs):
    keys = list(nrgs.keys())
    keys.sort()

    for key in keys:
        print("%25s : %12.8f" % (key, nrgs[key]))

system.setProperty("switchingFunction", HarmonicSwitchingFunction(10*angstrom, 9.5*angstrom))

printEnergies(system.energies())
 
print("\nEnergy with respect to cutoff length\n")
print("  Distance   Group    Shifted    ReactionField  ")

for i in range(10,501,5):
    x = i*0.1

    switchfunc = HarmonicSwitchingFunction(x*angstrom, (x-0.5)*angstrom)
    system.setProperty("switchingFunction", switchfunc)

    print("%12.8f  %12.8f  %12.8f  %12.8f" % (x, system.energy(group_coul).value(),
              system.energy(shift_coul).value(), system.energy(field_coul).value()))




