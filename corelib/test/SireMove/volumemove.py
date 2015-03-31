
from Sire.Mol import *
from Sire.Move import *
from Sire.FF import *
from Sire.MM import *
from Sire.System import *
from Sire.Units import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.IO import *
from Sire.Qt import *

timer = QTime()

# create a regular lattice of 125 LJ atoms (5 x 5 x 5)

###########
#
# This script has now been validated and shown to produce a
# trajectory that has the same average energy and volume
# as the trajectory modelled using ProtoMS. This validates
# that Sire is producing a correct NPT ensemble.
#
# The only downside is that ProtoMS produced 100M steps
# in 140 minutes, while Sire took 540 minutes...
#  (so Sire was nearly 4 times slower...)
#
#

mols = []

ljs = AtomLJs( [ LJParameter(2.0*angstrom, 2.0*kcal_per_mol) ] )

for i in range(0,5):
   for j in range(0,5):
       for k in range(0,5):
       
           editmol = EditMol("LJ")
           editmol.add( Atom("LJ", "H", Vector(2*i,2*j,2*k) ) )
           
           mol = Molecule(editmol)
           mol.setProperty("ljs", ljs)
           
           mols.append(mol)
           
PDB().write(mols, "test0000.pdb")

space = PeriodicBox( (0,0,0), (10,10,10) )

switchfunc = HarmonicSwitchingFunction( 4.5, 4.0 )

ljff = InterLJFF(space, switchfunc)
ljff.add(mols)

ffields = ForceFields()
ffields.add(ljff)

all_mols = MoleculeGroup("all", mols)

system = System(all_mols, ffields)
system.setSpace(space)

mc = RigidBodyMC( UniformSampler(all_mols) )
mc.setMaximumTranslation( 0.1 * angstrom )

volmc = VolumeMove( MapAsMolecules(all_mols) )
volmc.setVolumeChangingFunction( UniformVolumeChange(50 * angstrom3) )

moves = WeightedMoves()

moves.add(mc, 125)
moves.add(volmc, 1)

for i in range(1,1001):
    print("Running block %d" % i)
    timer.start()
    moves = system.run(moves, 100000)
    
    print("Took %d ms" % timer.elapsed())
    print("Energy = %f | Volume = %f" % (system.forceFields().energy(),    
                                         system.info().space().volume()))
                                         
    print("MC Accept %d Reject %d" % (moves.moves()[0].clone().nAccepted(),
                                      moves.moves()[0].clone().nRejected()))
    
    print("VMC Accept %d Reject %d" % (moves.moves()[1].clone().nAccepted(),
                                       moves.moves()[1].clone().nRejected()))
    
    
    
    PDB().write(system.forceFields().molecules(), "test%0004d.pdb" % i)
