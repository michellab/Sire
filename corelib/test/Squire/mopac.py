
from Sire.IO import *
from Sire.Squire import *
from Sire.Units import *

from Sire.Qt import *

water = PDB().readMolecule("test/io/sb2.pdb")

qmff = QMFF()

mopac = Mopac()

print(mopac.chargeCommandFile(water), "\n")

qmff.setQuantumProgram( mopac )

qmff.add(water)

print(qmff.energyCommandFile(), "\n")

t = QTime()

t.start()
nrg = qmff.energy()
ms = t.elapsed()

print("\nMopac energy = %s kcal mol-1" % nrg.to(kcal_per_mol))
print("Took %d ms" % ms)

t.start()
chgs = mopac.calculateCharges(water, {})
ms = t.elapsed()

print(chgs)
print(chgs.array())
print("Took %d ms" % ms)

