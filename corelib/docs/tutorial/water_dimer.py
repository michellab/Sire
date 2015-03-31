
# Import the Sire Molecule library. This library contains
# everything needed to represent and manipulate atoms
# and molecules
from Sire.Mol import *

# Import the Sire Input/Output library. This library contains
# everything needed to load and save things to and from files
from Sire.IO import *

# Import the Sire Maths library. This library contains mathematical
# objects, e.g. vectors, angles, planes etc.
from Sire.Maths import *

# Use the PDB object from Sire.IO to load a PDB file
# containing a water dimer
water_dimer = PDB().read("input/water_dimer.pdb")

# water_dimer is a container containing the two
# water molecules in the dimer. Lets get hold of
# each water molecule
first_water = water_dimer[MolIdx(0)].molecule()
print("The first water molecule is %s" % first_water)

# Now the second water molecule
second_water = water_dimer[MolIdx(1)].molecule()
print("The second water molecule is %s" % second_water)

# Let's print out all of the data about all of the atoms 
# of the first water
print("\nFirst water molecule:")
for i in range(0, first_water.nAtoms()):
    atom = first_water.atom(AtomIdx(i))
    print("%s: %s %s" % ( atom.name(),
                          atom.property("element"),
                          atom.property("coordinates") ))

# Let's do the same to the second water
print("\nSecond water molecule:")
for i in range(0, second_water.nAtoms()):
    atom = second_water.atom(AtomIdx(i))
    print("%s: %s %s" % ( atom.name(),
                          atom.property("element"),
                          atom.property("coordinates") ))

# Let's calculate the distance between the center of masses
# of the two water molecules
com_first_water = first_water.evaluate().centerOfMass()
com_second_water = second_water.evaluate().centerOfMass()

distance = Vector.distance(com_first_water, com_second_water)
print("\nThe distance between the centers of mass of the waters is %s A" % distance)
