from Sire.System import *
from Sire.Vol import *
from Sire.Base import *
from Sire.FF import *

sys = System()

from Sire.Maths import *

box0 = PeriodicBox( Vector(10.0,10.0,10.0) )
box1 = PeriodicBox( Vector(20.0,20.0,20.0) )

print(box0)
print(box0.volume())
print(box1.volume())

from Sire.MM import *
sys.add( InterCLJFF("cljff") )

print(sys)
print(sys.property("space"))
print(sys.userProperties().propertyKeys())
print(sys.builtinProperties().propertyKeys())

sys.setProperty( "space0", LinkToProperty("space", FFIdx(0)) )

print(sys.property("space0"))

sys.setProperty("space0", box0)

print(sys.property("space"))

sys.setProperty("space1", box1)

sys.setProperty("combined_space", CombineSpaces("space0", "space1"))

print(sys.properties().propertyKeys())

print(sys.property("combined_space"))
print(sys.property("combined_space").volume())

sys.setProperty("space0", PeriodicBox( Vector(5,5,5) ) )

print(sys.property("combined_space"))
print(sys.property("combined_space").volume())

sys.removeProperty("space0")

print(sys.properties().propertyKeys())

