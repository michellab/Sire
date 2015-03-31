from Sire.Base import *
from Sire.Mol import *
from Sire.IO import *

#load a TIP3P water molecule...
tip3p = PDB().read("test/geometry/tip3p.pdb")[0]

#add the interatomic bonds within this molecule...
tip3p.addAutoBonds()

#add the 'M' atom (initially at the origin)
tip3p[0].add("M03",Vector(0.0))

for atm in tip3p.atoms():
    print(atm)

#we need to set the coordinates of the M atom...

#first set the O-M bond length to 0.3 A
res = tip3p[0]
res.set(res.bond("M03","O00"), 0.3)

#now set the M-O-H angle to half the value of the H-O-H angle
hoh = res.angle("H01","O00","H02")
halfhoh = 0.5 * hoh.size()

res.set(res.angle("H01","O00","M03"), halfhoh, AbsFromMass())

#now set the M-H-O-H dihedral to 0
res.set( res.dihedral("M03","H01","O00","H02"), Angle(0.0) )

PDB().write(tip3p,"test.pdb")

AI = AtomIndex
hoh = tip3p.angle(AI("H01",1),AI("O00",1),AI("H02",1))
hom = tip3p.angle(AI("H01",1),AI("O00",1),AI("M03",1))

print(hoh)
print(hom)

for atm in tip3p.atoms():
    print(atm)


print("Done!")
