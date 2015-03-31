
from Sire.Base import *
from Sire.MM import *
from Sire.IO import *
from Sire.FF import *
from Sire.Units import *
from Sire.Maths import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Qt import *

t = QTime()

cljff = InterCLJFF()
fast_cljff = FastInterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

fast_cljff.setSwitchingFunction(switchfunc)
fast_cljff.setPatching( BoxPatching(vol, 15*angstrom) )

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
fast_cljff.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    cljff.add(mol) 
    fast_cljff.add(mol)   

solvent = MoleculeGroup("solvent", cljff.molecules())

delta = vol.dimensions()
print("Space == %s" % delta)

#for i in range(-1,1):
#    for j in range(-1,1):
#        for k in range(-1,1):
#            if i == 0 and j == 0 and k == 0:
#                continue
#
#            print "Adding molecules to box %d,%d,%d..." % (i,j,k)
#
#            add_mols = Molecules()
#
#            for imol in range(0,solvent.nMolecules()):
#                mol = solvent.moleculeAt(imol).molecule()
#                mol = mol.edit().renumber().commit()
#                mol = mol.move().translate( Vector( i*delta.x(), j*delta.y(), k*delta.z() ) ).commit()
#                add_mols.add(mol)
#
#            cljff.add(add_mols)
#            fast_cljff.add(add_mols)

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

mols = cljff.molecules()

#for molnum in mols.molNums():
#    mol = mols[molnum].molecule()
#    charges = mol.property("charge").array()
#    ljs = mol.property("LJ").array()
#
#    print charges, ljs

t.start()
cljff.packCoordinates()
ms = t.elapsed()
print("Packing the coordinates took %d ms" % ms)

#get the benchmark times
benchmark = 0.000001 * FlopsMark.benchmark() + 0.0001
benchmark_sum = 0.000001 * FlopsMark.benchmarkSum() + 0.0001
benchmark_quot = 0.000001 * FlopsMark.benchmarkQuotient() + 0.0001
benchmark_prod = 0.000001 * FlopsMark.benchmarkProduct() + 0.0001

print("\nThis machine can run at; %.1f MFLOPS for sum, %.1f MFLOPS for sum+product," % \
                (benchmark_sum, benchmark_prod))
print("%.1f MFLOPS for sum+quotient and %.1f MFLOPS for sum+product+sqrt.\n" % \
                (benchmark_quot, benchmark))

print(cljff.property("space"))

for i in range(0,1):
    t.start()
    cljff.mustNowRecalculateFromScratch()
    before_energy = FlopsMark()
    nrg = cljff.energy()
    after_energy = FlopsMark()
    ms = t.elapsed()
    print(nrg)
    print(nrg.value())

    mflops = 0.000001 * (after_energy - before_energy)

    print("Took %d ms. " % ms, end=' ')
    print("Speed is at least %.1f MFLOPS" % mflops)

    for j in range(0,before_energy.nThreads()):
        mflops_j = 0.000001 * (after_energy[j] - before_energy[j])
        print("%.1f MFLOPS for thread %d " % (mflops_j, j), end=' ')

    print("\n", end=' ')

    print("(This is %.2f %% of the benchmark  (%.2f %%, %.2f %%, %.2f %%))" % \
             ( 100 * (mflops / benchmark), 100 * (mflops / benchmark_quot), \
               100 * (mflops / benchmark_sum), 100 * (mflops / benchmark_prod) ))

    print(fast_cljff.patching())
    t.start()
    fast_cljff.mustNowRecalculateFromScratch()
    before_energy = FlopsMark()
    nrg = fast_cljff.energy()
    after_energy = FlopsMark()
    ms = t.elapsed()
    print(nrg)
    print(nrg.value())

    mflops = 0.000001 * (after_energy - before_energy)

    print("FastInterCLJFF Took %d ms. " % ms, end=' ')
    print("Speed is at least %.1f MFLOPS" % mflops)

    for j in range(0,before_energy.nThreads()):
        mflops_j = 0.000001 * (after_energy[j] - before_energy[j])
        print("%.1f MFLOPS for thread %d " % (mflops_j, j), end=' ')

    print("\n", end=' ')
    

print("Done!")

#now how long does it take to calculate a change in energy?
mols = cljff.molecules()
molnums = mols.molNums()

mol0 = mols[molnums[0]]

print(mol0.evaluate().center())

newmol = mol0.move().translate( Vector(1,0,0) ).commit()

print(mol0.evaluate().center())

cljff.update( newmol )
fast_cljff.update( newmol )

t.start()
nrg = cljff.energy()
ms1 = t.elapsed()

t.start()
nrg = fast_cljff.energy()
ms2 = t.elapsed()

print(ms1, ms2)

print("SLOW %s  (%f)" % (cljff.energy(), cljff.energy().value()))
print("FAST %s  (%f)" % (fast_cljff.energy(), fast_cljff.energy().value()))

cljff.update( mol0 )
fast_cljff.update( mol0 )

print("SLOW %s  (%f)" % (cljff.energy(), cljff.energy().value()))
print("FAST %s  (%f)" % (fast_cljff.energy(), fast_cljff.energy().value()))

print("Move two...")

mol1 = mols[molnums[1]]

print(mol1.evaluate().center())

mol1 = mol1.move().translate( Vector(1,0,0) ).commit()

print(mol1.evaluate().center())

mol0 = mol0.move().translate( Vector(2,0,0) ).commit()

print(mol0.evaluate().center())

cljff.update(mol0)
cljff.update(mol1)

fast_cljff.update(mol0)
fast_cljff.update(mol1)

t.start()
nrg1 = cljff.energy()
ms1 = t.elapsed()

t.start()
nrg2 = fast_cljff.energy()
ms2 = t.elapsed()

print("REAL %f kcal mol-1, VERSUS %f jcal mol-1" \
           % (nrg1.to(kcal_per_mol), nrg2.to(kcal_per_mol)))

print("TOOK %d versus %d milliseconds" % (ms1, ms2))
      
fast_cljff.mustNowRecalculateFromScratch();
print(fast_cljff.energy())

