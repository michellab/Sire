from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.FF import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *
from Sire.System import *
from Sire.Move import *
from Sire.Cluster import *
from Sire.Stream import *

import sys

t = QTime()

cljff = InterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("test/io/water.pdb")

print("Read in %d molecules!" % mols.nMolecules())

i = 0

t.start()
mol = mols.moleculeAt(0).molecule()

mol = mol.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.1550*kcal_per_mol)).molecule() \
                .atom( AtomName("H01") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("H02") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("M03") ) \
                    .setProperty("charge", -1.04 * mod_electron).molecule() \
         .commit()

charges = mol.property("charge")
ljs = mol.property("LJ")

cljff.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    cljff.add(mol)

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

system = System()

system.add(cljff)

print("Initial energy = %s" % system.energy())

data = save(system)

system = load(data)

print("Saved energy = %s" % system.energy())

mc = RigidBodyMC(cljff.group(MGIdx(0)))

moves = SameMoves(mc)

#give a lambda coordinate that turns off the coulomb energy
lam = Symbol("lambda")

system.setComponent( system.totalComponent(), cljff.components().lj() + \
                                              lam * cljff.components().coulomb() )

system.setComponent( lam, 0.0 )
print(system.energies())

system.setComponent( lam, 1.0 )
print(system.energies())

#create 5 replicas that map from lambda=0 to lambda=1
replicas = Replicas(system, 5)
replicas.setSubMoves(moves)
replicas.setNSubMoves(1000)

replicas.setLambdaComponent(lam)

for i in range(0, 5):
    replicas.setLambdaValue( i, i*0.0025 )

for i in range(0,5):
    print(i, replicas[i].lambdaValue(), replicas[i].energy(), \
             replicas[i].subSystem().constant(lam))

data = save(replicas)

replicas = load(data)

for i in range(0,5):
    print(i, replicas[i].lambdaValue(), replicas[i].energy())

repmove = RepExMove()

print("Running the replica exchange moves")
sim = SupraSim.run( replicas, repmove, 1, True )

replicas = sim.system()
repmove = sim.moves()

print("Accepted ",repmove[0].nAccepted(), "Rejected ",repmove[0].nRejected())

for i in range(0,5):
    print(i, replicas[i].lambdaValue(), replicas[i].energy(), \
             replicas[i].subSystem().constant(lam))
