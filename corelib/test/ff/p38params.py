
from SireMol import *
from SireIO import *
from SireFF import *

p38 = PDB().read("test/io/p38.pdb")[0]

ffdb = MMDB()

ProtoMS().read("test/ff/opls96.ff", ffdb)
ProtoMS().read("test/ff/opls96-residues.ff", ffdb)
ProtoMS().read("test/ff/solvents.ff", ffdb)

from SirePy import *
t = QTime()

t.start()

p38params = ffdb.getMMAtomParameters(p38, LamPosition.LAM_ANY)

print("Total charge on molecule = %f" % p38params.totalCharge())

for residue in p38:
    print("Total charge on residue %s (%d) = %f" % (residue.name(), residue.number(), p38params.totalCharge(residue.number())))

ms = t.elapsed()

print("Parametisation took %d ms" % ms)
