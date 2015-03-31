
from Sire.System import *
from Sire.Mol import *
from Sire.IO import *
from Sire.MM import *
from Sire.FF import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Units import *

cljff = InterCLJFF("CLJFF")

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("test/io/water.pdb")
                                                
print("Read in %d molecules!" % mols.nMolecules())

i = 0

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
    mols.update(mol)

system = System()

system.add(cljff)

print("System energy equals...")
print(system.energy())

group0 = MoleculeGroup("group0")
group1 = MoleculeGroup("group1")

group0.add( mols.moleculeAt(100) )
group0.add( mols.moleculeAt(101) )

group1.add( mols.moleculeAt(102) )
group1.add( mols.moleculeAt(103) )
group1.add( mols.moleculeAt(104) )

cljff2 = InterGroupCLJFF("group_energy")
cljff2.add( group0, MGIdx(0) )
cljff2.add( group1, MGIdx(1) )
cljff2.setSpace(vol)
cljff2.setSwitchingFunction(switchfunc)

print(cljff2.energy())

system.add(cljff2)
system.add(group0)
system.add(group1)

print(system.energies())

nrgmon = EnergyMonitor(group0, group1)

print(nrgmon)

nrgmon.monitor(system)

cnrgs = nrgmon.coulombEnergies()
ljnrgs = nrgmon.ljEnergies()

for i in range(0,cnrgs.nRows()):
    for j in range(0,cnrgs.nColumns()):
        print(i,j,cnrgs(i,j).average(),ljnrgs(i,j).average())

print(nrgmon.views0())
print(nrgmon.views1())

nrgmon.monitor(system)

cnrgs = nrgmon.coulombEnergies()
ljnrgs = nrgmon.ljEnergies()

for i in range(0,cnrgs.nRows()):  
    for j in range(0,cnrgs.nColumns()):
        print(i,j,cnrgs(i,j).average(),ljnrgs(i,j).average())
        print(i,j,cnrgs(i,j).nSamples(),ljnrgs(i,j).nSamples())

