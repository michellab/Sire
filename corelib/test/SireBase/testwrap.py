
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Base import *

water = PDB().readMolecule("test/io/water.pdb")

center = water.evaluate().center()

dblarray = [ 1.0,2,3,4,5 ]
intarray = [ 1,2,3,4,5 ]
vecarray = [ Vector(1), Vector(2), Vector(3) ]
strarray = [ "cat", "dog", "fish" ]


water = water.edit().setProperty("center", wrap(center)) \
                    .setProperty("dblarray", wrap(dblarray)) \
                    .setProperty("intarray", wrap(intarray)) \
                    .setProperty("vecarray", wrap(vecarray)) \
                    .setProperty("strarray", wrap(strarray)) \
                    .setProperty("type", wrap("ligand")) \
                    .setProperty("alpha", wrap(0.5)) \
                    .setProperty("copies", wrap(1)).commit()


print((water.property("center")))
print((water.property("dblarray")))
print((water.property("intarray")))
print((water.property("vecarray")))
print((water.property("strarray")))
print((water.property("type")))
print((water.property("alpha")))
print((water.property("copies")))

print((water.properties()))

