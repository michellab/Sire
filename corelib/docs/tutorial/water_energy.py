
from Sire.Maths import *

# The Sire.MM library contains everything needed to calculate
# molecular mechanics (MM) energies
from Sire.MM import *

# Import the Sire.Stream library so that we can load
# back the water dimer from "water_dimer.s3"
import Sire.Stream

# Load the parameterised water dimer from the Sire Saved Stream (s3) file
(first_water, second_water) = Sire.Stream.load("water_dimer.s3")

# Create an intermolecular coulomb and Lennard Jones forcefield (InterCLJFF).
# This is a container that calculates the coulomb and Lennard Jones energies
# between all contained molecules
cljff = InterCLJFF("clj")

# Add both water molecules to the forcefield
cljff.add( first_water )
cljff.add( second_water )

# The total (coulomb plus LJ) energy is returned using the "energy()" function
total_nrg = cljff.energy()

print("The total interaction energy between the dimer is %s" % total_nrg)

# You can get the coulomb and LJ components using;
coul_nrg = cljff.energy( cljff.components().coulomb() )
lj_nrg = cljff.energy( cljff.components().lj() )

print("The coulomb energy is %s. The LJ energy is %s." % (coul_nrg, lj_nrg))

# Let's calculate the energy as a function of the intermolecular
# distance between the water dimer. First, get the vector between
# the center of masses of the waters
com_vector = second_water.evaluate().centerOfMass() - \
             first_water.evaluate().centerOfMass()

# How far apart are they now?
com_distance = com_vector.length()
print("\nThe water molecules are separated by %s A" % com_distance)

# Now translate them so that they are separated by 2 A
com_vector = com_vector.normalise()
delta_vector = (2 - com_distance) * com_vector

print("\nTo move the water molecules to be separated by 2 A, we " \
      "need to translate the second water by %s" % delta_vector)

# Actually translate the second water molecule. This uses the 
# "moves().translate()" function
second_water = second_water.move().translate(delta_vector).commit()

# Check that this has worked...
com_distance = Vector.distance( first_water.evaluate().centerOfMass(),
                                second_water.evaluate().centerOfMass() )

print("The updated distance is indeed %s A" % com_distance)

# Update the forcefield with the new version of "second_water"
cljff.update(second_water)

# Print all energy components using the "energies" function
print(cljff.energies())

# Now loop over every distance from 2 - 4 A, and calculate the energies
print("\nDistance   Energies")
for i in range(0,10):
    second_water = second_water.move().translate( 0.2*com_vector ).commit()
    cljff.update(second_water)
    print("%f  %s" % ( Vector.distance( first_water.evaluate().centerOfMass(), \
                                        second_water.evaluate().centerOfMass() ),
                       cljff.energies() ))

