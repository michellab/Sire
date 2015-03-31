from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.ID import *

protein = PDB().readMolecule("test/io/cox2.pdb")

mass = 0             
com = Vector(0)

for atom in protein.atoms( AtomName("CA", CaseInsensitive) ):
    com += atom.property("coordinates") * atom.property("element").mass().value()
    mass += atom.property("element").mass().value()

print("Center of mass = %s" % (com / mass))

