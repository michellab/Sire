
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
from Sire.Stream import *

import sys

t = QTime()

cljff = InterCLJFF("cljff")

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

solvent = MoleculeGroup("solvent")

mols = PDB().read("test/io/water.pdb")
                                                
print("Read in %d molecules!" % mols.nMolecules())

i = 0

t.start()
mol = mols.moleculeAt(0).molecule()

mol = mol.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.15500*kcal_per_mol)).molecule() \
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

solvent.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    solvent.add(mol)

    cljff.add(mol)

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

system = System()
system.add(solvent)
system.add(cljff)

t.start()
print("Initial energy = %s" % system.energy())
print("(took %d ms)" % t.elapsed())
print(system.property("space"))
print(system.property("switchingFunction"))

print(system.groupNumbers())
print(system.groupNames())
print(system.energies())

mc = RigidBodyMC(solvent)

moves = SameMoves(mc)

moves.setGenerator( RanGenerator(42) )

print("Running 10000 moves without saving the trajectory...")
t.start()
system = moves.move(system, 10000, True)
ms = t.elapsed()
print("Done! (took %d ms)" % ms)

system.add( "trajectory", TrajectoryMonitor(solvent), 1000 )

print("Running 10000 moves with saving the trajectory...")
t.start()
system = moves.move(system, 10000, True)
ms = t.elapsed()
print("Done! (took %d ms)" % ms)

print("Writing the trajectory to disk...")
t.start()

system[ MonitorName("trajectory") ].writeToDisk("tempXXXXXX.pdb")

ms = t.elapsed()
print("Took %d ms" % ms)

print("Final energy = %s" % system.energy())

system.mustNowRecalculateFromScratch();

print("Are we sure? = %s" % system.energy())

mc = moves.moves()[0]

print("nAccepted() == %d, nRejected() == %d  (%f %%)" % (mc.nAccepted(), \
                            mc.nRejected(), 100 * mc.acceptanceRatio()))
