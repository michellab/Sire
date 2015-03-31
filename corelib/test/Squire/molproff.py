
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.FF import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Base import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *
from Sire.Squire import *
from Sire.Cluster import *
from Sire.System import *
from Sire.Move import *

import os
import sys

timer = QTime()

#read in all of the molecules
print("Loading and parameterising the molecules...")
timer.start()
mols = PDB().read("test/io/water.pdb")

mol = mols.moleculeAt(0).molecule()

mol = mol.edit().atom( AtomName("H01") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("H02") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("M03") ) \
                    .setProperty("charge", -1.04 * mod_electron).molecule() \
                .atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter( 3.15365 * angstrom, 
                                                    0.1550 * kcal_per_mol ) ) \
         .molecule().commit()

charges = mol.property("charge")
ljs = mol.property("LJ")

ms = timer.elapsed()
print("... took %d ms" % ms)

#specify the space in which the molecules are placed
space = Cartesian()

space = PeriodicBox(Vector(-18.3854,-18.66855,-18.4445), \
                    Vector( 18.3854, 18.66855, 18.4445))

#specify the type of switching function to use
switchfunc = HarmonicSwitchingFunction(80.0*angstrom)
switchfunc = HarmonicSwitchingFunction(15.0*angstrom, 14.5*angstrom)

#create a forcefield for the molecules
qmff = QMFF("QMFF")

qmff.setProperty("space", space)

molpro = Molpro()
molpro.setEnvironment("HELP", "THIS IS SOME HELP")

qmff.setProperty("quantum program", molpro)

qmff.add(mols.moleculeAt(0))

print("QM energy in current thread")

qmnrg = qmff.energy()

print(qmnrg)

qmmmff = QMMMFF("qmmmff")

qmmmff.setSpace(space)
qmmmff.setSwitchingFunction(switchfunc)
qmmmff.setQuantumProgram(molpro)

qmmmff.add(mols.moleculeAt(0), MGIdx(0))

mmff = InterGroupCLJFF("mmff")
mmff.setSpace(space)
mmff.setSwitchingFunction(switchfunc)

mmff.add(mol, MGIdx(0))

for i in range(1, 5): #mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    qmmmff.add(mol, MGIdx(1))
    mmff.add(mol, MGIdx(1))

qmmmnrg = qmmmff.energy()
mmnrg = mmff.energy( mmff.components().coulomb() )

FILE = open("molpro.cmd", "w")

print(qmmmff.energyCommandFile(), file=FILE)

FILE.close()

sys.exit(0)

print(qmmmnrg)
print(qmmmnrg - qmnrg)
print(mmnrg)
print(qmmmnrg - qmnrg - mmnrg)
print(mmff.energy( mmff.components().lj() ))


for i in range(0,100):
     print("Step %d" % i)
     qmmmff.mustNowRecalculateFromScratch()
     print(qmmmff.energy())

sys.exit(0)

system = System()

system.add(qmmmff)

print("Initial energy = %s" % system.energy())

system.mustNowRecalculateFromScratch()

mc = RigidBodyMC(qmmmff.group(MGIdx(0)))

moves = SameMoves(mc)

mtsmc = MTSMC()

print("Running 5 moves using MPI")
nodes = MPINodes()
node = nodes.getFreeNode()
print("node rank = %d of %d" % (node.rank(), nodes.nNodes()))

sim = Simulation.run(node, system, moves, 5)

print("Job submitted. Waiting...")

sim.wait()

print("Job complete!")

system = sim.system()

print("Final energy = %s" % system.energy())

system.mustNowRecalculateFromScratch();

print("Are we sure? = %s" % system.energy())

mc = sim.moves().moves()[0]

print("nAccepted() == %d, nRejected() == %d  (%f %%)" % (mc.nAccepted(), \
                            mc.nRejected(), 100 * mc.acceptanceRatio()))

print("Took %d ms" % ms)
