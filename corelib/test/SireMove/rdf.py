
from Sire.Mol import *
from Sire.Move import *
from Sire.FF import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.CAS import *
from Sire.IO import *
from Sire.Qt import *

timer = QTime()

solvent = PDB().read("test/io/water.pdb")

solvent_group = MoleculeGroup("solvent",solvent)

solute = PDB().read("test/io/tip4p.pdb")[0]

rdf = RDF(0,15,45)

oxygen = CGAtomID(CutGroupID(0), AtomID(0))
hydrogen0 = CGAtomID(CutGroupID(0), AtomID(1))
hydrogen1 = CGAtomID(CutGroupID(0), AtomID(2))

for mol in solvent:
    rdf += Vector.distance( solute[oxygen], 
                            mol[oxygen] )
    
for point in rdf.normalise():
    print("%f  %f" % point)


