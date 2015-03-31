
from Sire.Base import *
from Sire.MM import *
from Sire.IO import *
from Sire.FF import *
from Sire.System import *
from Sire.Move import *
from Sire.Units import *
from Sire.Maths import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Qt import *

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

t.start()
nrg = system.energy()
ms = t.elapsed()

print("Energy = %f kcal mol-1 - took %d ms" % (nrg.to(kcal_per_mol), ms))

closemols1 = CloseMols( Vector(0,0,0), cljff[MGIdx(0)], 10 )
closemols2 = CloseMols( Vector(0,0,0), cljff[MGIdx(0)], 20 )
closemols3 = CloseMols( Vector(0,0,0), cljff[MGIdx(0)], 1 )

closemols1.update(system)
closemols2.update(system)
closemols3.update(system)

mols1 = list(closemols1.closeMolecules().keys())
mols1.sort()

mols2 = list(closemols2.closeMolecules().keys())
mols2.sort()

mols3 = list(closemols3.closeMolecules().keys())
mols3.sort()

print("CLOSEMOLS1: ",mols1)
print("CLOSEMOLS2: ",mols2)
print("CLOSEMOLS3: ",mols3)

move = RigidBodyMC( cljff[MGIdx(0)] )

move.move(system, 100)

print("New energy = %f kcal mol-1" % (system.energy().to(kcal_per_mol)))
