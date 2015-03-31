
from Sire.Base import *

p = PackedArray2D_double_( [ [1,2,3], [4,5,6,7] ] )

print(p)
print((p[0]))
print((p[1]))
print((p[0][0]))

p.append( [10,11,12,13,14,15] )

print(p)
print((p[2]))

p.remove(0)

print(p)

p = PackedArray2D_double_(p,p)

print(p)

p.remove(3)
p.remove(1)

print(p)

