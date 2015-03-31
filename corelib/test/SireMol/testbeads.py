
from Sire.IO import *
from Sire.Mol import *

print("Loading protein...")
mol = PDB().readMolecule("test/io/p38.pdb")

print("Creating one bead per residue...")
beads = Beads(mol)

print("Printing out the beads...")
print(beads)

