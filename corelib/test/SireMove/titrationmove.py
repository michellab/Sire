
from Sire.IO import *
from Sire.Mol import *
from Sire.Move import *
from Sire.System import *
from Sire.MM import *
from Sire.Units import *

import sys

(molecules, space) = Amber().readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")

water = molecules[MolNum(1)].molecule()

water = water.edit() \
             .atom(AtomName("O")) \
             .setProperty("PDB-atom-name", "O") \
             .setProperty("element", Element("O")) \
             .molecule().atom(AtomName("H1")) \
             .setProperty("PDB-atom-name", "H1") \
             .setProperty("element", Element("H")) \
             .molecule().atom(AtomName("H2")) \
             .setProperty("PDB-atom-name", "H2") \
             .setProperty("element", Element("H")) \
             .molecule().commit()

sodium = molecules[MolNum(1)].molecule()

sodium = sodium.edit() \
               .atom(AtomName("O")) \
               .setProperty("PDB-atom-name", "Na") \
               .setProperty("element", Element("Na")) \
               .setProperty("charge", 1*mod_electron) \
               .setProperty("LJ", LJParameter(3.3284*angstrom, 0.0030*kcal_per_mol)) \
               .molecule().atom(AtomName("H1")) \
               .setProperty("PDB-atom-name", "Xx") \
               .setProperty("charge", 0*mod_electron) \
               .setProperty("LJ", LJParameter.dummy()) \
               .setProperty("element", Element(0)) \
               .molecule().atom(AtomName("H2")) \
               .setProperty("charge", 0*mod_electron) \
               .setProperty("LJ", LJParameter.dummy()) \
               .setProperty("PDB-atom-name", "Xx") \
               .setProperty("element", Element(0)) \
               .molecule().commit()

chloride = molecules[MolNum(1)].molecule()

chloride = chloride.edit() \
                   .atom(AtomName("O")) \
                   .setProperty("PDB-atom-name", "Cl") \
                   .setProperty("charge", -1*mod_electron) \
                   .setProperty("LJ", LJParameter(4.4010*angstrom, 0.1000*kcal_per_mol)) \
                   .setProperty("element", Element("Cl")) \
                   .molecule().atom(AtomName("H1")) \
                   .setProperty("PDB-atom-name", "Xx") \
                   .setProperty("charge", 0*mod_electron) \
                   .setProperty("LJ", LJParameter.dummy()) \
                   .setProperty("element", Element(0)) \
                   .molecule().atom(AtomName("H2")) \
                   .setProperty("PDB-atom-name", "Xx") \
                   .setProperty("charge", 0*mod_electron) \
                   .setProperty("LJ", LJParameter.dummy()) \
                   .setProperty("element", Element(0)) \
                   .molecule().commit()

titrator = Titrator()
titrator.setPositiveTemplate(sodium)
titrator.setNegativeTemplate(chloride)
titrator.setNeutralTemplate(water)

system = System()

cljff = InterCLJFF("cljff")
cljff.add(molecules)

solvent = MoleculeGroup("solvent")
solvent.add(molecules)

titrator.setMoleculeGroup(solvent)

system.add(cljff)
system.add(solvent)
system.setProperty("space", space)

print("Initialising the ions...")
titrator.applyTo(system)
print("Randomising the location of ions...")
titrator.randomiseCharge(3)
titrator.applyTo(system)
print("System is ready for simulation :-)")

PDB().write(system.molecules(), "test0000.pdb")

move = TitrationMove()
move.setTemperature( 25*celsius )

moves = WeightedMoves()
moves.add(move, 1)

move = RigidBodyMC(solvent)
moves.add(move, 1)

print("Start: %s" % system.energies())

for i in range(1,11):
    system = moves.move(system, 1000, False)
    print("%5d: %s" % (i, system.energies()))
    PDB().write(system.molecules(), "test%0004d.pdb" % i)

print(moves)

