
from Sire.Mol import *
from Sire.IO import *

mols = PDB().read("test/geometry/tip3p.pdb")

tip3p = mols[0]

mols = PDB().read("test/geometry/tip4p.pdb")
tip4p = mols[0]
    
print(tip3p)
print(tip4p)

print("OLD TIP3P")
for atom in tip3p[0]:
    print(atom)

tip3p[0].applyTemplate(tip4p[0],TmplType(TmplTypeEnum.MATCHATOMS))

print("TIP4P")
for atom in tip4p[0]:
    print(atom)
    
print("NEW TIP3P")
for atom in tip3p[0]:
    print(atom)

print("Done!")
