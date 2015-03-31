
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.Qt import *
from Sire.Units import *
from Sire.Move import *
from Sire.Maths import *

import sys

t = QTime()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)

mols = PDB().read("test/io/water.pdb")
                                                
print("Read in %d molecules!" % mols.nMolecules())

PDB().write(mols, "test00.pdb")

for i in range(0,5):
    new_vol = vol.setVolume( 1.05 * vol.volume() )

    print("update...")
    for j in range(0, mols.nMolecules()):
        mols.update( mols[MolIdx(j)].move().changeSpace(vol, new_vol) )

    print("write...")
    PDB().write(mols, "test0%d.pdb" % (i+1))
    print("DONE")

    vol = new_vol

