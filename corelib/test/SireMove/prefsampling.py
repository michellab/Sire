
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
from Sire.Cluster import *
from Sire.Move import *

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
mol0 = mols.moleculeAt(0).molecule()

mol0 = mol0.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.1550*kcal_per_mol)).molecule() \
                .atom( AtomName("H01") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("H02") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("M03") ) \
                    .setProperty("charge", -1.04 * mod_electron).molecule() \
         .commit()

charges = mol0.property("charge")
ljs = mol0.property("LJ")

cljff.add(mol0)

mol1 = mols.moleculeAt(1).molecule()

mol1 = mol1.edit().rename("T4P") \
                  .setProperty("charge", charges) \
                  .setProperty("LJ", ljs).commit()

cljff.add(mol1)

mol2 = mols.moleculeAt(2).molecule()

mol2 = mol2.edit().rename("T4P") \
                  .setProperty("charge", charges) \
                  .setProperty("LJ", ljs).commit()

cljff.add(mol2)

system = System()

system.add(cljff)

print("Initial energy = %s" % system.energy())

sampler = PrefSampler(mol0, cljff[MGIdx(0)], 200*angstrom2)
sampler.updateFrom(system)

mol0 = PartialMolecule(mol0)
mol1 = PartialMolecule(mol1)
mol2 = PartialMolecule(mol2)

p0 = sampler.probabilityOf(mol0)
p1 = sampler.probabilityOf(mol1)
p2 = sampler.probabilityOf(mol2)

print(p0, p1, p2, p0+p1+p2) 

mol1 = mol1.move().translate( Vector(1,0,0) ).commit()

system.update(mol1)
sampler.updateFrom(system)

p0 = sampler.probabilityOf(mol0)
p1 = sampler.probabilityOf(mol1)
p2 = sampler.probabilityOf(mol2)

print(p0, p1, p2, p0+p1+p2)

mol0 = mol0.move().translate( Vector(1,0,0) ).commit()

system.update(mol0)
sampler.updateFrom(system)

p0 = sampler.probabilityOf(mol0)
p1 = sampler.probabilityOf(mol1)
p2 = sampler.probabilityOf(mol2)
                                                   
print(p0, p1, p2, p0+p1+p2)

