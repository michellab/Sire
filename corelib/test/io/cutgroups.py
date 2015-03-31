from SireMol import *
from SireIO import *

mollist = PDB().read("test/geometry/tip4p.pdb")

if (mollist.size() == 0):
    print("Error - we have not loaded a molecule!")
    
tip4p = mollist[0]
tip4p.addAutoBonds()
tip4p.setName("TIP4P")

#create a cutgroup
cg = Molecule.createFrom(tip4p, MolType.RIGIDMOL)
print(cg.toString())

for atom in cg.atoms():
    print(atom.toString())
    
print("Goodbye!")
