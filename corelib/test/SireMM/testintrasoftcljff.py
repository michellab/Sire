
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Units import *

print("Loading the molecules...")
(molecules, space) = Amber().readCrdTop("test/io/SYSTEM.crd", "test/io/SYSTEM.top")

# get the protein, which is the first molecule
protein = molecules[MolNum(2)].molecule()

print("The protein has %d atoms." % protein.nAtoms())
print("The protein had %d residues." % protein.nResidues())

print("Building the forcefields...")

# test with different coul and LJ cutoffs
ffs = [ IntraCLJFF("a"), IntraCLJFF("b"),
        IntraCLJFF("c"), IntraCLJFF("d") ]

soft_ffs = [ IntraSoftCLJFF("a"), IntraSoftCLJFF("b"), 
             IntraSoftCLJFF("c"), IntraSoftCLJFF("d") ]

short = 15 * angstrom
long = 25 * angstrom
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
print("Testing IntraCLJFF and IntraSoftCLJFF")
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
print("Testing IntraGroupCLJFF and IntraGroupSoftCLJFF")
print("===============================================")

# test with different coul and LJ cutoffs
ffs = [ IntraGroupCLJFF("a"), IntraGroupCLJFF("b"),
        IntraGroupCLJFF("c"), IntraGroupCLJFF("d") ]
    
soft_ffs = [ IntraGroupSoftCLJFF("a"), IntraGroupSoftCLJFF("b"),
             IntraGroupSoftCLJFF("c"), IntraGroupSoftCLJFF("d") ]

residues = protein.residues()

for i in range(0,int(residues.count() / 2)):
    for j in range(0,4):
        ffs[j].add(residues[i], MGIdx(0))
        soft_ffs[j].add(residues[i], MGIdx(0))

for i in range(int(residues.count() / 2), residues.count()):
    for j in range(0,4):
        ffs[j].add(residues[i], MGIdx(1))
        soft_ffs[j].add(residues[i], MGIdx(1))

for i in range(0,4):
    ffs[i].setSwitchingFunction(switchfuncs[i])
    soft_ffs[i].setSwitchingFunction(switchfuncs[i])
            
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

