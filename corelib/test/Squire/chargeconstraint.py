
from Sire.IO import *
from Sire.Squire import *
from Sire.System import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Maths import *
from Sire.Move import *
from Sire.Vol import *
from Sire.Units import *

import os

mol = PDB().readMolecule("test/io/methanol.pdb")

mol = mol.edit().rename("methanol").commit()

protoms = ProtoMS("%s/Work/ProtoMS/protoms2" % os.getenv("HOME"))

protoms.addParameterFile("test/ff/methanol.par")
protoms.addParameterFile("test/ff/solvents.ff")

mol = protoms.parameterise(mol, ProtoMS.SOLUTE)

print(mol.property("charge").array())

qmff = QMFF("MopacFF")
qmff.setQuantumProgram( Mopac() )

qmff.add(mol)

intraff = InternalFF("IntraFF")
intraff.add(mol)

intraclj = IntraCLJFF("IntraCLJ")
intraclj.add(mol)

print(intraff.energy())
print(intraclj.energy())

print(intraff.energy() + intraclj.energy())

print(qmff.energy())

solute = MoleculeGroup("solute")
solute.add(mol)

system = System()

system.add(qmff)
system.add(intraff)
system.add(intraclj)

system.add(solute)

print(system.energy())

chg_constraint = QMChargeConstraint(solute)
chg_constraint.setChargeCalculator( AM1BCC() )

system.add(chg_constraint)

print(system.constraintsSatisfied())

print(system.energy())

print(system.constraintsSatisfied())

mol = system[ MGIdx(0) ][ mol.number() ].molecule()

print(mol.property("charge").array())

water = PDB().read("test/io/water.pdb")

tip4p = water[ MolIdx(0) ].molecule()
tip4p = protoms.parameterise(tip4p, ProtoMS.SOLVENT)

chgs = tip4p.property("charge")
ljs = tip4p.property("LJ")

for i in range(0, water.nMolecules()):
    tip4p = water[ MolIdx(i) ].molecule()
    tip4p = tip4p.edit().setProperty("charge", chgs) \
                        .setProperty("LJ", ljs) \
                 .commit()

    water.update(tip4p)

words = open("test/io/water.xsc", "r").readlines()[0].split()

space = PeriodicBox( Vector( float(words[0]), float(words[1]), float(words[2]) ),
                     Vector( float(words[3]), float(words[4]), float(words[5]) ) )

cljff = InterGroupCLJFF("solute-solvent")

cljff.add( mol, MGIdx(0) )
cljff.add( water, MGIdx(1) )

system.add(cljff)

waterff = InterCLJFF("solvent")
waterff.add( water )

system.add(waterff)

all = MoleculeGroup("all")
all.add(mol)
all.add(water)

system.add(all)

system.setProperty("space", space)

rbmc = RigidBodyMC(all)
zmatmc = ZMatMove(solute)

moves = WeightedMoves()
moves.add(rbmc, 500)
moves.add(zmatmc)

for i in range(0,100):
    print("MOVE")
    system = moves.move( system, 500, True )
    print("ENERGY = %f kcal mol-1" % (system.energy().to(kcal_per_mol)))
    print(moves)

    mol = system[ solute.name() ][ MolIdx(0) ].molecule()
    print("CHARGES ",mol.property("charge"))

