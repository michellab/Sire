
from Sire.Vol import *
from Sire.IO import *

import copy

coords = PDB().readMolecule("test/io/sb2.pdb") \
              .molecule().property("coordinates").array()

coords2 = copy.copy(coords)

print((coords == coords2))

print(coords.nCoordGroups(), coords.nCoords())

coords.remove(1,2)

assert( coords.nCoordGroups() == coords2.nCoordGroups()-2 )
assert( coords.nCoords() == coords2.nCoords() - 
                            coords2[1].count() - coords2[2].count() )

print(coords.nCoordGroups(), coords.nCoords())
print(coords2.nCoordGroups(), coords2.nCoords())

print((coords == coords2))

ngroups = coords.nCoordGroups()
ncoords = coords.nCoords()

coords.append(coords2)

assert( coords.nCoordGroups() == ngroups + coords2.nCoordGroups() )
assert( coords.nCoords() == ncoords + coords2.nCoords() )

print(coords.nCoordGroups(), coords.nCoords())

while coords.nCoordGroups() != coords2.nCoordGroups():
    coords.remove(0)

assert( coords == coords2 )

print(coords.nCoordGroups(), coords.nCoords())

coords.append(coords)
coords.append(coords)

assert( coords.nCoordGroups() == coords2.nCoordGroups() * 4 )
assert( coords.nCoords() == coords2.nCoords() * 4 )

for i in range(0,coords2.nCoordGroups()):
    assert( coords[i] == coords[i+coords2.nCoordGroups()] )
    assert( coords[i] == coords[i+2*coords2.nCoordGroups()] )
    assert( coords[i] == coords[i+3*coords2.nCoordGroups()] )

print(coords.nCoordGroups(), coords.nCoords())

coords.remove(0, coords.nCoordGroups())

print(coords.nCoordGroups(), coords.nCoords())

assert( coords.isEmpty() )

