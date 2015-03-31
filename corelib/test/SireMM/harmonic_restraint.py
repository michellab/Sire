
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.Units import *
from Sire.Maths import *

water = PDB().readMolecule("test/io/water.pdb")

atom = water.atom( AtomName("O00") )
point = atom.property("coordinates")

k = 25 * kcal_per_mol / (angstrom*angstrom)

# construct a restraint that is a harmonic potential
# between 'point' and 'atom', using the passed force
# constant
restraint = DistanceRestraint.harmonic(atom, point, k)

# the above line is a quick short-hand for the general
# DistanceRestraint constructor, which allows any 
# energy function to be used, based on the distance
# between the two points, represented by the symbol
# "r" - e.g. we could create the same restraint using;

# r = DistanceRestraint.r()
# restraint = DistanceRestraint(atom, point, k.value() * r**2) 

# (note that DistanceRestraint.r() is just a convenient
# shorthand for 'r = Symbol("r")' - it is there so that
# you don't have to remember which symbol is used to 
# represent the distance between the two points of the
# restraint

print(restraint.energy())

water = water.move().translate( Vector(1,0,0) ).commit()

restraint.update(water)

print(restraint.energy())
