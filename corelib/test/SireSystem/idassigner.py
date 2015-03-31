
from Sire.System import *
from Sire.IO import *
from Sire.FF import *
from Sire.Maths import *
from Sire.Mol import *
from Sire.Vol import *

mols = PDB().read("test/io/water.pdb")

solvent = MoleculeGroup("solvent", mols)

system = System()
system.add(solvent)
system.setProperty("space", Cartesian())

points = [Vector(10,10,10), solvent.moleculeAt(0).atom(AtomName("O00")), 
          Vector(0,0,0), Vector(5,5,5)]

print(solvent.moleculeAt(0).molecule().number())

idassigner = IDAssigner(points, solvent)

points = idassigner.points()

idassigner.update(system)

mols = idassigner.identifiedMolecules()

for i in range(0,len(points)):
    print(points[i], points[i].point(), mols[i].number(), mols[i].evaluate().center())

