
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Units import *

print("Loading the molecules...")
(molecules, space) = Amber().readCrdTop("test/io/SYSTEM.crd", "test/io/SYSTEM.top")

# get the protein, which is the first molecule
protein = molecules[MolNum(2)].molecule()

water = Molecules()
for i in range(3, molecules.nMolecules()):
    water.add( molecules[MolNum(i)].molecule() )

print("The protein has %d atoms." % protein.nAtoms())
print("There solvent has %d molecules." % water.nMolecules())

print("Building the forcefields...")

# test with different coul and LJ cutoffs
ffs = [ InterCLJFF("a"), InterCLJFF("b"),
        InterCLJFF("c"), InterCLJFF("d") ]

soft_ffs = [ InterSoftCLJFF("a"), InterSoftCLJFF("b"), 
             InterSoftCLJFF("c"), InterSoftCLJFF("d") ]

short = 10 * angstrom
long = 20 * angstrom
feather = 0.5 * angstrom
rf_dielectric = 78.3

switchfuncs = [ HarmonicSwitchingFunction(long, long-feather, long, long-feather),
                HarmonicSwitchingFunction(short, short-feather, short, short-feather),
                HarmonicSwitchingFunction(long, long-feather, short, short-feather),
                HarmonicSwitchingFunction(short, short-feather, long, long-feather) ]

for i in range(0,4):
    ffs[i].setSwitchingFunction(switchfuncs[i])
    soft_ffs[i].setSwitchingFunction(switchfuncs[i])

    ffs[i].add(protein)
    soft_ffs[i].add(protein)

    ffs[i].add(water)
    soft_ffs[i].add(water)

    ffs[i].setProperty("space", space)
    soft_ffs[i].setProperty("space", space)

    ffs[i].setUseReactionField(True)
    ffs[i].setReactionFieldDielectric(rf_dielectric)

    soft_ffs[i].setUseReactionField(True)
    soft_ffs[i].setReactionFieldDielectric(rf_dielectric)

    soft_ffs[i].setShiftDelta(2.0)
    soft_ffs[i].setCoulombPower(0)

def printEnergies(ffs, soft_ffs):
    print("Hard[L:L  S:S  L:S  S:L] Soft[L:L  S:S  L:S  S:L]")

    for i in range(0,11):
        alpha = 0.1 * i

        for ff in soft_ffs:
            ff.setAlpha(alpha)

        print("COUL: %8.5f " % alpha, end=' ')

        for j in range(0,4):
            print("%9.3f " % ffs[j].energy( ffs[j].components().coulomb() ).to(kcal_per_mol), end=' ')

        for j in range(0,4):
            print("%9.3f " % soft_ffs[j].energy( soft_ffs[j].components().coulomb() ).to(kcal_per_mol), end=' ')

        print("\nLJ  : %8.5f " % alpha, end=' ')

        for j in range(0,4):
            print("%9.3f " % ffs[j].energy( ffs[j].components().lj() ).to(kcal_per_mol), end=' ')
                
        for j in range(0,4):
            print("%9.3f " % soft_ffs[j].energy( soft_ffs[j].components().lj() ).to(kcal_per_mol), end=' ')

        print("\n", end=' ')

print("=====================================")
print("Testing InterCLJFF and InterSoftCLJFF")
print("=====================================")

print("Reaction field cutoff")
printEnergies(ffs, soft_ffs)

for i in range(0,4):
    ffs[i].setUseGroupCutoff(True)
    soft_ffs[i].setUseGroupCutoff(True)

print("\nGroup feathered cutoff")
printEnergies(ffs, soft_ffs)

for i in range(0,4):
    ffs[i].setShiftElectrostatics(True)
    soft_ffs[i].setShiftElectrostatics(True)

print("\nElectrostatic shifted cutoff")
printEnergies(ffs, soft_ffs)

for i in range(0,4):
    ffs[i].setUseAtomisticCutoff(True)
    soft_ffs[i].setUseAtomisticCutoff(True)

print("\nAtomistic cutoff")
printEnergies(ffs, soft_ffs)

print("\n===============================================")
print("Testing InterGroupCLJFF and InterGroupSoftCLJFF")
print("===============================================")

# test with different coul and LJ cutoffs
ffs = [ InterGroupCLJFF("a"), InterGroupCLJFF("b"),
        InterGroupCLJFF("c"), InterGroupCLJFF("d") ]
    
soft_ffs = [ InterGroupSoftCLJFF("a"), InterGroupSoftCLJFF("b"),
             InterGroupSoftCLJFF("c"), InterGroupSoftCLJFF("d") ]

for i in range(0,4):
    ffs[i].setSwitchingFunction(switchfuncs[i])
    soft_ffs[i].setSwitchingFunction(switchfuncs[i])

    ffs[i].add(protein, MGIdx(0))
    soft_ffs[i].add(protein, MGIdx(0))
    ffs[i].add(water, MGIdx(1))
    soft_ffs[i].add(water, MGIdx(1))
            
    ffs[i].setUseReactionField(True)
    ffs[i].setReactionFieldDielectric(rf_dielectric)

    soft_ffs[i].setUseReactionField(True)
    soft_ffs[i].setReactionFieldDielectric(rf_dielectric)

    soft_ffs[i].setShiftDelta(2.0)
    soft_ffs[i].setCoulombPower(0)

print("Reaction field cutoff")
printEnergies(ffs, soft_ffs)
    
for i in range(0,4):
    ffs[i].setUseGroupCutoff(True)
    soft_ffs[i].setUseGroupCutoff(True)

print("\nGroup feathered cutoff")
printEnergies(ffs, soft_ffs)
            
for i in range(0,4):
    ffs[i].setShiftElectrostatics(True)
    soft_ffs[i].setShiftElectrostatics(True)
    
print("\nElectrostatic shifted cutoff")   
printEnergies(ffs, soft_ffs)
    
for i in range(0,4):
    ffs[i].setUseAtomisticCutoff(True)
    soft_ffs[i].setUseAtomisticCutoff(True)
        
print("\nAtomistic cutoff")
printEnergies(ffs, soft_ffs)

