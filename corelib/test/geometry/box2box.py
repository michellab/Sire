from Sire.Mol import *
from Sire.IO import *

tip4p = PDB().read("test/geometry/tip4p.pdb")[0]

watbox = PDB().read("test/geometry/tip3pbox.pdb")

for tip3p in watbox:
    tip3p[0].applyTemplate(tip4p[0],TmplType(TmplTypeEnum.MATCHATOMS))

PDB().write(watbox,"test.pdb")

print("Done!")
