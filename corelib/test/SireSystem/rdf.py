
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
from Sire.MPI import *

import sys

t = QTime()

cljff = InterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15, 14.5)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("test/io/water.pdb")
                                                
print("Read in %d molecules!" % mols.nMolecules())

i = 0

t.start()
mol = mols.moleculeAt(0).molecule()

mol = mol.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363,  \
                                                   0.1550)).molecule() \
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

    mol = mol.edit().rename( MolName("T4P") ) \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    cljff.add(mol)

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

system = System()

rdf = RDFMonitor( 0*angstrom, 10*angstrom, 0.05*angstrom )
rdf.add( MolIdx(0) + AtomName("O00"), AtomName("O00") )

rdf2 = RDFMonitor( 0*angstrom, 10*angstrom, 0.05*angstrom )
rdf2.add( MolIdx(0) + AtomName("O00"), AtomName("H01") )
rdf2.add( MolIdx(0) + AtomName("O00"), AtomName("H02") )

system.add(cljff)
system.add("O-O RDF", rdf)
system.add("O-H RDF", rdf2)

t.start()
system.collectStats()

ms = t.elapsed()

print("Collecting stats took %d ms" % ms)

rdf = system.monitor( MonitorName("O-O RDF") )

print("\nO-O RDF")

for i in range(0,rdf.nBins()):
    print(rdf[i].middle().to(angstrom), rdf[i].value())

rdf = system.monitor( MonitorName("O-H RDF") )

print("\nO-H RDF")

for i in range(0,rdf.nBins()):
    print(rdf[i].middle().to(angstrom), rdf[i].value())

