
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

import Sire.Stream

t = QTime()

cljff = InterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("test/io/water.pdb")
                                                
print(("Read in %d molecules!" % mols.nMolecules()))

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
print(("Parameterised all of the water molecules (in %d ms)!" % ms))

system = System()

system.add(cljff)

lam = Symbol("lambda")

system.setComponent( lam, 0.2 )
system.setComponent( system.totalComponent(), lam * cljff.components().total() )

mc = RigidBodyMC(cljff.group(MGIdx(0)))

moves = SameMoves(mc)

def testStream(c):
    t.start()

    data = Sire.Stream.save(c)

    ms = t.elapsed()

    print(("Streaming %s took %d ms" % (c.what(), ms)))
    print(("%s takes up %d bytes" % (c.what(),data.size())))

    t.start()

    header = Sire.Stream.getDataHeader(data)
    print((header.dataType()))
    print((header.requiredLibraries()))

    print((header.createdBy()))
    print((header.createdWhere()))

    print((header.requiredMemory()))
    print((header.compressionRatio()))
    print((header.digest()))
    print((header.repository()))
    print((header.buildVersion()))
    print((header.systemInfo()))

    c2 = Sire.Stream.load(data)

    ms = t.elapsed()
  
    print(("Reading the data took %d ms" % ms))
    print(c)
    print(c2)

testStream(system)

data = Sire.Stream.save(system)

print("Probing the system...")
print((system.energy()))
print((system.energies()))

system = Sire.Stream.load(data)

print((system.energy()))
print((system.energies()))

print("\nGetting data info...")

t.start()
Sire.Stream.save( system, "test/SireStream/tmp_testdata.sire" )
ms = t.elapsed()

print(("Saving a system to a file took %d ms" % ms))

t.start()
system = Sire.Stream.load( "test/SireStream/tmp_testdata.sire" )
ms = t.elapsed()

print(("Reading a system from a file took %d ms" % ms))

header = Sire.Stream.getDataHeader( "test/SireStream/tmp_testdata.sire" )
print((header.dataType()))
print((header.requiredLibraries()))

print((header.createdBy()))
print((header.createdWhere()))

print((header.requiredMemory()))
print((header.compressionRatio()))
print((header.digest()))
print((header.repository()))
print((header.buildVersion()))
print((header.systemInfo()))

print((system.energies()))

system.mustNowRecalculateFromScratch()

print((system.energies()))

