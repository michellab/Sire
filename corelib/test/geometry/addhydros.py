from Sire.Mol import *
from Sire.IO import *

mols = PDB().read("test/geometry/dioxin-no-hyd.pdb")

if (mols.count() > 0):
    dioxin_noh = mols[0]

mols = PDB().read("test/geometry/dioxin.pdb")
if (mols.count() > 0):
    dioxin = mols[0]
    
print("DIOXIN, no hydrogens")
for atom in dioxin_noh[0]:
    print(atom)

dioxin_noh[0].applyTemplate(dioxin[0],TmplType(TmplTypeEnum.MATCHATOMS))

print("DIOXIN")
for atom in dioxin[0]:
    print(atom)
    
print("NEW DIOXIN, added hydrogens!")
for atom in dioxin_noh[0]:
    print(atom)

PDB().write(dioxin_noh,"test.pdb")

print("Done!")

