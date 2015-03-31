
from Sire.Mol import *
from Sire.IO import *
from Sire.Maths import *

tip4p = PDB().read("test/geometry/tip4p.pdb")[0]

print("Self same")
print(tip4p.signature())
assert(tip4p.signature() == tip4p.signature())

tip2 = PDB().read("test/geometry/tip4p.pdb")[0]

print("Copy same")
assert(tip2 == tip4p)
assert(tip2.signature() == tip4p.signature())

tip2.translate( Vector(1.0,0.0,0.0) )

print("Sig same")
assert( not (tip2 == tip4p) )
assert(tip2.signature() == tip4p.signature())

tip2.remove( AtomIndex("O00",1) )

print("Sig different")
assert( not (tip2 == tip4p) )
assert( not (tip2.signature() == tip4p.signature()) )

tip3 = EditMol(tip4p.name())
tip3.addResidue(1,"WTR")
tip3.add( Atom("H01",1) )
tip3.add( Atom("H02",1) )
tip3.add( Atom("M03",1) )
tip3.add( Atom("O00",1) )

print("Manual same?")
assert( not (tip3 == tip4p) )
assert( tip3.signature() == tip4p.signature() )


print("Done!")
