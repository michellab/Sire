
from Sire.MM import *
from Sire.Mol import *
from Sire.IO import *
from Sire.Maths import *
from Sire.Units import *

# The Sire.Vol library contains everything needed to define
# spaces (volumes), e.g. periodic boundary condition boxes
from Sire.Vol import *  

# The Sire.Move library contains all of the different types
# of move (Monte Carlo and molecular dynamics)
from Sire.Move import *

# The Sire.System library contains everything needed to define
# a simulation system
from Sire.System import *

# The Sire.CAS library contains a complete Computer Algebra System
# that lets you build and manipulate algebraic expressions
from Sire.CAS import *

import Sire.Stream

# First, lets load up the parameterised water dimer from the s3 file
(first_water, second_water) = Sire.Stream.load("water_dimer.s3")

# Now we create the InterCLJFF and add the two water molecules
cljff = InterCLJFF("cljff")
cljff.add(first_water)
cljff.add(second_water)

# Sire performs simulations on simulation System objects.
# A System groups together all of the molecules, forcefields etc.
# into a single unit
system = System()
system.add(cljff)

# Here we create the equation used to calculate the total energy of the
# system. We define a symbol, lambda, which we use here to scale up and
# down the coulomb component of cljff
lam = Symbol("lambda")
total_nrg = cljff.components().lj() + lam * cljff.components().coulomb()

# Here we tell the system to use our equation to calculate the total
# energy, and we set the value of lambda to 0.5
system.setComponent( system.totalComponent(), total_nrg )
system.setConstant( lam, 0.5 )

# Now we create a MoleculeGroup that groups together all of the 
# molecules to be moved
mobile_mols = MoleculeGroup("mobile_molecules")
mobile_mols.add(first_water)
mobile_mols.add(second_water)

# We add the molecule group to the system
system.add(mobile_mols)

# We create a periodic boundaries space, and give it to the system
space = PeriodicBox( Vector(5,5,5) )
system.setProperty( "space", space )

# We also add a SpaceWrapper that ensures any moves on mobile_mols
# will wrap the molecules back into the simulation space (centered
# on the origin)
system.add( SpaceWrapper(Vector(0), mobile_mols) )

# Now we create the move object that performs rigid body translation
# and rotation Monte Carlo moves on the water dimer (molecules in 
# the mobile_mols MoleculeGroup)
move = RigidBodyMC(mobile_mols)

# We set the move temperature, and also the maximum amounts by
# which we can translate or rotate molecules. Note that we have
# to specify units (e.g. celsius, angstrom, degrees. It is perfectly
# ok to use fahrenheit, kelvin, nanometers, radians etc. if you want)
move.setTemperature( 25*celsius )
move.setMaximumTranslation(0.5 * angstrom)
move.setMaximumRotation(5 * degrees)

# We now group all of the moves to be performed into a single Moves
# object. In this case, a WeightedMoves will draw moves randomly
# according to their weight
moves = WeightedMoves()
moves.add( move, 1 )

# Lets perform 100 moves. The moves are performed on a copy of 'system',
# with the updated version of 'system' after the moves returned by this
# function
print("Running 100 moves...")
new_system = moves.move(system, 100, True)

# Now lets run a simulation, writing out a PDB trajectory
PDB().write(system.molecules(), "output000.pdb")
print(system.energies())

# Here we run 10 blocks of 1000 moves, printing out the 
# energies and a PDB of the coordinates after each block
for i in range(1,11):
    system = moves.move(system, 1000, True)
    print("%d: %s" % (i, system.energies()))
    PDB().write(system.molecules(), "output%003d.pdb" % i)

# Finally, we print out information about how many moves
# were accepted and rejected.
print("Move information:")
print(moves)

