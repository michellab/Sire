
from Sire.Tools import WSRC
from Sire.IO import *
from Sire.Maths import *
from Sire.Mol import *

import Sire.Stream

import os

(waters, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

print(space.dimensions())

ligand = Sire.Stream.load("../io/ligand.s3")
ligand = ligand.move().translate( -ligand.evaluate().center() ).commit()
ligand = ligand.move().translate( 0.5 * space.dimensions() ).commit()

overlaps = WSRC.getOverlapWaters(ligand, waters)

print(overlaps)

swap = MoleculeGroup("swap")

for overlap in overlaps:
    swap.add( overlap.molecule() )

PDB().write(ligand, "ligand.pdb")
PDB().write(swap, "waters.pdb")

command = "rm ligand.pdb waters.pdb"
os.system(command)

if __name__ == "__main__":
    print("OK")
