
from Sire.CAS import *

f = Function("f")
g = Function("g")

x = Symbol("x")
y = Symbol("y")

ex = f(x) + g(x)

print(ex)

ex = f( g(x) )

print(ex)

print(ex.substitute( g(x) == pow(x,2) ))
print(ex.substitute( f(g(x)) == pow(x,2) ))

ex = ex.differentiate(x)

print(ex)

print(ex.substitute( f(g(x)) == pow(x,2) ))

print("Done!")
