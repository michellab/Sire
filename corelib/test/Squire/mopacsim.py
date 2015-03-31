
from Sire.IO import *
from Sire.Squire import *
from Sire.Units import *
from Sire.System import *
from Sire.Move import *
from Sire.Mol import *

from Sire.Qt import *

import os

sb2 = PDB().readMolecule("test/io/sb2.pdb")

sb2 = sb2.edit().rename("SB2").commit()

protoms_dir = os.getenv("HOME") + "/Work/ProtoMS"

protoms = ProtoMS( "%s/protoms2" % protoms_dir )

protoms.addParameterFile("%s/parameter/amber99.ff" % protoms_dir )
protoms.addParameterFile("%s/parameter/gaff.ff" % protoms_dir )
protoms.addParameterFile("test/ff/sb2.ff")

sb2 = protoms.parameterise(sb2, ProtoMS.SOLUTE)

qmff = QMFF("MopacFF")
mopac = Mopac()
qmff.setQuantumProgram( mopac )

qmff.add(sb2)

system = System()

system.add(qmff)

zmat_move = ZMatMove(qmff[MGIdx(0)])

moves = SameMoves(zmat_move)

import Sire.Stream
Sire.Stream.save( (system,moves), "test/Squire/mopacsim.s3" )

print("Energy before == %f kcal mol-1" % system.energy().to(kcal_per_mol))
system = moves.move(system)
print("Energy after == %f kcal mol-1" % system.energy().to(kcal_per_mol))
