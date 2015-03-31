
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

cljff_a = InterCLJFF("cljff_a")
cljff_b = InterCLJFF("cljff_b")
cljff_a_b = InterGroupCLJFF("cljff_a_b")

ffs = [ cljff, cljff_a, cljff_b, cljff_a_b ]

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

for ff in ffs:
    ff.setSpace(vol)
    ff.setSwitchingFunction(switchfunc)

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
cljff_a.add(mol)
cljff_a_b.add(mol, MGIdx(0))

solvent.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
                    .setProperty("center", wrap(mol.evaluate().center())) \
             .commit()

    solvent.add(mol)

    cljff.add(mol)

    if (i % 2 == 0):
        cljff_a.add(mol)
        cljff_a_b.add( mol, MGIdx(0) )
    else:
        cljff_b.add(mol)
        cljff_a_b.add( mol, MGIdx(1) )

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

system2 = System()
system2.add(solvent)

system2.add(cljff_a)
system2.add(cljff_b)
system2.add(cljff_a_b)

print(system.groupNumbers())
print(system.groupNames())
print(system.energies())
print(system2.groupNumbers())
print(system2.groupNames())
print(system2.energies())

data = save(system2)

t.start()
print("Other initial energy = %s" % system2.energy())
print("(took %d ms)" % t.elapsed())
print(system2.property("space"))
print(system2.property("switchingFunction"))

mc = RigidBodyMC(solvent)

moves = SameMoves(mc)

t.start()
moves.setGenerator( RanGenerator(42) )

for i in range(0,10):
    system = moves.move(system, 1, False)
    print("Energy = %s" % system.energy())

ms = t.elapsed()

print("Done!")

print("Final energy = %s" % system.energy())

system.mustNowRecalculateFromScratch();

print("Are we sure? = %s" % system.energy())

mc = moves.moves()[0]

print("nAccepted() == %d, nRejected() == %d  (%f %%)" % (mc.nAccepted(), \
                            mc.nRejected(), 100 * mc.acceptanceRatio()))

moves = SameMoves(mc)
moves.clearStatistics()

t.start()
moves.setGenerator( RanGenerator(42) )

print("Running 5000 moves")
system2 = moves.move(system2, 5000, False)
print("Energy = %s" % system2.energy())

ms = t.elapsed()

print("Done! - took %d ms" % ms)

print("Final energy (2) = %s" % system2.energy())

system2.mustNowRecalculateFromScratch();

print("Are we sure? (2) = %s" % system2.energy())

mc = moves.moves()[0]

print("nAccepted() == %d, nRejected() == %d  (%f %%)" % (mc.nAccepted(), \
                            mc.nRejected(), 100 * mc.acceptanceRatio()))

system3 = load(data)
print(system3.energies())

moves = SameMoves(mc)
moves.clearStatistics()

t.start()
moves.setGenerator( RanGenerator(42) ) 

for i in range(0,10):
    system3 = moves.move(system3, 1, False)
    print("Energy = %s" % system3.energy())

ms = t.elapsed()
                            
print("Done!")

print("Final energy = %s" % system3.energy())

system3.mustNowRecalculateFromScratch();

print("Are we sure? = %s" % system3.energy())

mc = moves.moves()[0]

print("nAccepted() == %d, nRejected() == %d  (%f %%)" % (mc.nAccepted(), \
                            mc.nRejected(), 100 * mc.acceptanceRatio()))

# Check that there is no drift in the center of the molecules
for molnum in system2.molNums():
    water = system2[molnum].molecule()

    if water.hasProperty("center"):
        center = water.property("center")
        eval_center = water.evaluate().center()
        distance = Vector.distance(center, eval_center)

        if distance > 0.1:
            print("WARNING: Drift in atom center: %s %s %s" % (center, eval_center, distance))

