
import Sire.Stream

from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Base import *
from Sire.Maths import *
from Sire.Units import *

ligand = Sire.Stream.load("test/io/osel.s3")
ligand = ligand.edit().setProperty("center", wrap(ligand.evaluate().center())).commit()

print("Original center = %s" % ligand.property("center"))

intraff = InternalFF("intraff")
intraff.add(ligand)

intraclj = IntraCLJFF("intraclj")
intraclj.add(ligand)

system = System()
system.add(intraff)
system.add(intraclj)

mols = MoleculeGroup("mols")
mols.add(ligand)
system.add(mols)

intramove = InternalMove(mols)
rbmove = RigidBodyMC(mols)
rbmove.setMaximumTranslation(0*angstrom)

moves = WeightedMoves()
moves.add(intramove, 1)
moves.add(rbmove, 1)

for i in range(0,10):
    system = moves.move(system, 250, False)
    print("Completed 250 moves...")
    PDB().write(system.molecules(), "test%0004d.pdb" % (i+1))

    ligand = system[ligand.number()].molecule()

    print("New center = %s" % ligand.property("center"))

print("Complete")
