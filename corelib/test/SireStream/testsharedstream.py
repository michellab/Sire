from Sire.Base import *
from Sire.Stream import *
from Sire.Qt import *

props = Properties()

p = StringProperty("mieow")

props.setProperty( "cat", p )
props.setMetadata( "cat", p )
props.setProperty( "tiger", p )
props.setMetadata( "tiger", p )

data = save(props)

print(("Properties takes up %d bytes" % data.size()))

print("Loading the properties...")
p2 = load(data)
print("...complete!")

print(props)
print(p2)

print((props.property("cat")))
print((p2.property("cat")))

