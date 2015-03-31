
from Sire.System import *
from Sire.Mol import *
from Sire.IO import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Vol import *

lam = Symbol("lambda")
alpha = Symbol("alpha")

d_lambda = Delta( lam, 0.0, 1.0 )
d_alpha = Delta( alpha, 0.5, 0.6 )

print(d_lambda)
print(d_alpha)
print((d_lambda + d_alpha))

d_lambda2 = Delta( lam, 1.0, 0.0 )
d_alpha2 = Delta( alpha, 0.6, 0.5 )

print(d_lambda2)
print(d_alpha2)

try:
    print((d_lambda + d_lambda2))
    assert( False )
except:
    print("Cannot combine %s and %s (expected)" % (d_lambda, d_lambda2))

try:
    print((d_alpha + d_alpha2))
    assert( False )
except:
    print("Cannot combine %s and %s (expected)" % (d_alpha, d_alpha2))

assert( (d_lambda + d_alpha).involves(lam) )
assert( (d_lambda + d_alpha).involves(alpha) )
assert( not d_lambda.involves(alpha) )
assert( not d_alpha.involves(lam) )

water = PDB().readMolecule("test/io/water.pdb")

d_water = Delta( water, water.move().translate( Vector(1,0,0) ).commit() )

print(d_water)
assert( d_water.involves(water) )

d_water2 = Delta( water, water.move().translate( Vector(2,0,0) ).commit() )
assert( d_water2.involves(water) )

print(d_water2)
print((d_water + d_water2))
print((d_water2 + d_water))
assert( (d_water + d_water2).involves(water) )

water2 = water.edit().renumber().commit()

d_water3 = Delta( water2, water2.move().translate( Vector(1,0,0) ).commit() )

print(d_water3)

assert( d_water3.involves(water2) )
assert( not d_water3.involves(water) )

print((d_water + d_water3))
big_delta = (d_water + d_lambda + d_alpha + d_water3)
print(big_delta)

assert( big_delta.involves(lam) )
assert( big_delta.involves(alpha) )
assert( big_delta.involves(water) )
assert( big_delta.involves(water2) )

print(Delta(water, Molecule()))
print(Delta(Molecule(), water))

assert( Delta(water, Molecule()).involves(water) )
assert( Delta(Molecule(), water).involves(water) )

print(Delta(Molecules(water), Molecules()))
print(Delta(Molecules(), Molecules(water)))

assert( Delta(Molecules(water), Molecules()).involves(water) )
assert( Delta(Molecules(), Molecules(water)).involves(water) )

print((Delta(water,Molecule()) + Delta(Molecule(),water)))
assert( not (Delta(water,Molecule()) + Delta(Molecule(),water)).involves(water) )

print(Delta(water, water2))
assert( Delta(water, water2).involves(water) )
assert( Delta(water, water2).involves(water2) )

print(Delta(water2, water))

null_delta = (Delta(water, water2) + Delta(water2, water))
print(null_delta)
assert( not null_delta.involves(water) )
assert( not null_delta.involves(water2) )
assert( null_delta.isEmpty() )
assert( not null_delta.isNull() )

water3 = water.move().translate( Vector(1,1,1) ).commit()
d_water = Delta(water, water3)
d_water2 = Delta(water3, water)
print(d_water)
print(d_water2)

try:
    print((d_water + d_water2))
    assert(False)
except:
    print("Cannot combine %s and %s (expected)" % (d_water, d_water2))

print(Delta(water, water))
assert( Delta(water, water).isEmpty() )

print((d_water + Delta(water, water)))
assert( (d_water + Delta(water, water)).isEmpty() )

try:
    print((d_water + Delta(water3, water3)))
    assert(False)
except:
    print("Cannot combine %s and %s (expected)" % (d_water, Delta(water3,water3)))

cart = Cartesian()
box = PeriodicBox(Vector(10,10,10))
idcon = IdentityConstraint()
pcon = PropertyConstraint()

space_delta = Delta("space", cart, box)
print(space_delta)
assert( space_delta.involves("space") )
assert( not space_delta.involves("constraint") )

con_delta = Delta("constraint", idcon, pcon)
print(con_delta)
assert( con_delta.involves("constraint") )
assert( not con_delta.involves("space") )

print((space_delta + con_delta))
assert( (space_delta+con_delta).involves("space") )
assert( (space_delta+con_delta).involves("constraint") )

space_delta1 = Delta("space", box, cart)
print(space_delta1)

con_delta1 = Delta("constraint", pcon, idcon)
print(con_delta1)

try:
    print((space_delta + space_delta1))
    assert(False)
except:
    print("Cannot combine %s and %s (expected)" % (space_delta, space_delta1))

try:
    print((con_delta1 + con_delta))
    assert(False)
except:
    print("Cannot combine %s and %s (expected)" % (con_delta1, con_delta))

print((big_delta + space_delta + con_delta))
