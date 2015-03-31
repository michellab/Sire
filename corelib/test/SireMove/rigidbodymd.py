
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

mols = PDB().read("test/io/water.pdb")
                                                
print("Read in %d molecules!" % mols.nMolecules())

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

cljff = InterCLJFF("water-water")

cljff.add(mol)

solvent = MoleculeGroup("solvent")
solvent.add(mol)

for i in range(1,7):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()

    solvent.add(mol)
    cljff.add(mol)

system = System()
system.add(solvent)
system.add(cljff)

print(system.energy())

rbmove = MolecularDynamics( solvent, DLMRigidBody(), 1*femtosecond )

#rbmove.setEnergyComponent( cljff.components().coulomb() )

PDB().write(system.molecules(), "test0000.pdb")

for i in range(1,1000):
    rbmove.move(system, 10)
    print(i, system.energy())
    print(rbmove.kineticEnergy(), (system.energy() + rbmove.kineticEnergy()))

    PDB().write(system.molecules(), "test%0004d.pdb" % i)

