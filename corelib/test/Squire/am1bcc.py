
from Sire.IO import *
from Sire.Squire import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Maths import *
from Sire.Units import *
from Sire.Qt import *

import os

mol = PDB().readMolecule("test/io/methanol.pdb")

mol = mol.edit().rename("methanol").commit()

protoms = ProtoMS("%s/Work/ProtoMS/protoms2" % os.getenv("HOME"))

protoms.addParameterFile("test/ff/methanol.par")
protoms.addParameterFile("test/ff/solvents.ff")

mol = protoms.parameterise(mol, ProtoMS.SOLUTE)

am1bcc = AM1BCC()

t = QTime()

t.start()
am1bcc_chgs = am1bcc(mol)
ms = t.elapsed()

print("Getting QM charges took %d ms" % ms)

print("MM charges == ",mol.property("charge").array())
print("QM charges == ",am1bcc_chgs.array())

water = PDB().read("test/io/water.pdb")

tip4p = water[ MolIdx(0) ].molecule()

tip4p = protoms.parameterise(tip4p, ProtoMS.SOLVENT)

chgs = tip4p.property("charge")
ljs = tip4p.property("LJ")

for i in range(0, water.nMolecules()):
    tip4p = water[ MolIdx(i) ].molecule()
    tip4p = tip4p.edit().setProperty("charge", chgs) \
                        .setProperty("LJ", ljs) \
                 .commit()

    water.update(tip4p)

cljff = InterGroupCLJFF()

cljff.add( mol, MGIdx(0) )
cljff.add( water, MGIdx(1) )

print("MM energy == %f kcal mol-1" % (cljff.energy().to(kcal_per_mol)))
print("  coulomb == %f, LJ == %f" % (cljff.energy(cljff.components().coulomb()).to(kcal_per_mol),
                                     cljff.energy(cljff.components().lj()).to(kcal_per_mol)))

newmol = mol.edit().setProperty("charge", am1bcc_chgs).commit()
cljff.update(newmol)

mol = mol.move().translate(Vector(10,0,0)).commit()

t.start()
print(am1bcc.mayChangeCharges(mol, newmol))
ms = t.elapsed()

print("Checking for need to recalculate charges took %d ms" % ms)

mol = mol.move().rotate( Quaternion(25*degrees, Vector(1,1,1)),  Vector(0,0,0) ).commit()

t.start()
print(am1bcc.mayChangeCharges(mol, newmol))
ms = t.elapsed()

print("Checking for need to recalculate charges took %d ms" % ms)

mol = mol.atom(AtomIdx(5)).move().translate( Vector(3,3,3) ).commit().molecule()

t.start()
print(am1bcc.mayChangeCharges(mol, newmol))
ms = t.elapsed()

print("Checking for need to recalculate charges took %d ms" % ms)

print("QM energy == %f kcal mol-1" % (cljff.energy().to(kcal_per_mol)))
print("  coulomb == %f, LJ == %f" % (cljff.energy(cljff.components().coulomb()).to(kcal_per_mol),
                                     cljff.energy(cljff.components().lj()).to(kcal_per_mol)))
