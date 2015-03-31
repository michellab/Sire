
from Sire.Vol import *
from Sire.Maths import *
from Sire.Units import *

r = RegularGrid( Vector(1,2,3), 5, 2 * angstrom )

print(r)
print(r.center())
print(r.gridSpacing())
print(r.points())

r = r.rotate( Quaternion( 32*degrees, Vector(1,0,0) ), r.center() )

print(r)
print(r.center())
print(r.gridSpacing())
print(r.points())

