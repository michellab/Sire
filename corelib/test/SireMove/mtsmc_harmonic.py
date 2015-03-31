
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
from Sire.Base import *

import sys

timer = QTime()

#read in the solute molecule
print("Loading the tip4p solute...")
tip4p = PDB().readMolecule("test/io/tip4p.pdb")

#specify the space in which the molecules are placed
space = PeriodicBox( (-18.3854,-18.66855,-18.4445),
                     ( 18.3854, 18.66855, 18.4445) )

hff = HarmonicFF(space)

tip4p.setProperty( "k", QVariant(10.0) )
tip4p.setProperty( "r0", QVariant.fromValue(tip4p.extract().geometricCenter()) )

hff.add( tip4p, { hff.parameters().k() : "k",
                  hff.parameters().r0() : "r0" } )

hff2 = HarmonicFF(space)

tip4p.setProperty( "k2", QVariant(7.5) )
hff2.add( tip4p, {hff2.parameters().k() : "k2",
                  hff2.parameters().r0() : "r0" } )

e_slow = FFExpression( "e_slow", hff.components().total() )
e_fast = FFExpression( "e_fast", hff2.components().total() )

solute = MoleculeGroup("solute")
solute.add(tip4p)

groups = MoleculeGroups(solute)
ffields = ForceFields()
ffields.add(hff)
ffields.add(hff2)

ffields.add(e_slow)
ffields.add(e_fast)

ffields.setTotal(e_slow)

monitors = SystemMonitors()
monitors.set( ffields.total(), MonitorEnergy(ffields.total()) )

system = System(groups, ffields, monitors)

mc = RigidBodyMC( UniformSampler(solute) )

mc.setTemperature( 25 * celsius )
mc.setMaximumTranslation( 0.3 * angstrom )

#for i in range(0,50):
#    moves = system.run(mc, 10000)
#
#    mc = moves.moves()[0].clone()
#    print "%d accepted, %d rejected, ratio = %f %%" % \
#             (mc.nAccepted(), mc.nRejected(), 100 * mc.acceptanceRatio())
#         
#    print system.monitors().monitor(ffields.total()).average()


mtsmc = MTSMC(mc, e_fast.function(), 500)
mtsmc.setEnergyComponent(e_slow.function())

for i in range(0,5000):
    moves = system.run(mtsmc, 20)

    mtsmc = moves.moves()[0].clone()
    print("%d accepted, %d rejected, ratio = %f %%" % \
             (mtsmc.nAccepted(), mtsmc.nRejected(), 100 * mtsmc.acceptanceRatio()))
         
    print("AVGENERGY: %f" % system.monitors().monitor(ffields.total()).average())
    print("ENERGY: %f" % system.forceFields().energy())
    sys.stdout.flush()
