
from Sire.IO import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *

amber = Amber()

(waters, space) = amber.readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")

group_cljff = InterCLJFF("group_cutoff")

shift_cljff = InterCLJFF("shift_cutoff")
shift_cljff.setShiftElectrostatics(True)

field_cljff = InterCLJFF("field_cutoff")
field_cljff.setUseReactionField(True)
field_cljff.setReactionFieldDielectric(78.0)

atom_cljff = InterCLJFF("atomistic_cutoff")
atom_cljff.setUseAtomisticCutoff(True)

group_coul = group_cljff.components().coulomb()
shift_coul = shift_cljff.components().coulomb()
field_coul = field_cljff.components().coulomb()
atom_coul = atom_cljff.components().coulomb()

forcefields = [ group_cljff, shift_cljff, field_cljff, atom_cljff ]

switchfunc = HarmonicSwitchingFunction(10*angstrom, 9.5*angstrom)

system = System()

for forcefield in forcefields:
    forcefield.add(waters)
    system.add(forcefield)

def printEnergies(nrgs):
    keys = list(nrgs.keys())
    keys.sort()

    for key in keys:
        print("%25s : %12.8f" % (key, nrgs[key]))

system.setProperty("space", space)
system.setProperty("switchingFunction", switchfunc)

printEnergies(system.energies())
 
print("\nEnergy with respect to cutoff length\n")
print("  Distance   Group    Shifted    ReactionField   Atomistic")

for i in range(10,200,5):
    x = i*0.1

    switchfunc = HarmonicSwitchingFunction(x*angstrom, (x-0.5)*angstrom)
    system.setProperty("switchingFunction", switchfunc)

    print("%12.8f  %12.8f  %12.8f  %12.8f  %12.8f" % (x, system.energy(group_coul).value(),
              system.energy(shift_coul).value(), system.energy(field_coul).value(),
              system.energy(atom_coul).value()))
