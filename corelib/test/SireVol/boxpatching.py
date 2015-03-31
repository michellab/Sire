
from Sire.Vol import *
from Sire.Maths import *

space = PeriodicBox( Vector(-15,-14,-13), Vector(15,14,13) )
patches = BoxPatching(space)

print(patches)

print(patches.patchIndex( Vector(-15,-14,-13) ))
print(patches.patchIndex( Vector(14.999,14,13) ))
print(patches.patchIndex( Vector(15,14,13) ))

