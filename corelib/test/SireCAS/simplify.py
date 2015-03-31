
from Sire.CAS import *

x = Symbol("x")
y = Symbol("y")
z = Symbol("z")

f = (1-z) * (x+y) + z * x

print(f)
print(f.expand(x))
print(f.simplify())

print(f( {x:1, y:2, z:3} ))
print(f.simplify()( {x:1, y:2, z:3} ))

g = (z*x) + (1-z) * (x+y)
print(g)
print(g.expand(x))
print(g.simplify())

print(g( {x:1, y:2, z:3} ))
print(g.simplify()( {x:1, y:2, z:3} ))

