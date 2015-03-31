
from Sire.Mol import *
import Sire.Maths

from Sire.IO import *

tip4p = EditMol("TIP4P")

tip4p.add( ResNum(1), "WTR" )

residue = tip4p[ ResNum(1) ]

print(residue)

residue.add("H01")
residue.add("H02")
residue.add("O00")
residue.add("M03")

print(residue)
print(tip4p)

#for atom in residue:
#    print atom.name()

residue.addBond("H01","O00")
residue.addBond("H02","O00")
residue.addBond("M03","O00")

residue.set(Bond("H01","O00",ResNum(1)), 0.96)
residue.set(Bond("H02","O00",ResNum(1)), 0.96)
residue.set(Bond("M03","O00",ResNum(1)), 0.3)

hoh = Sire.Maths.Angle.degrees(109.5)

residue.set(Angle("H01","O00","H02",ResNum(1)), hoh)

#check_angle = residue.measure("H01","O00","H02")
#print "Set angle to %f degrees, but is actually %f degrees." % (hoh.toDegrees(), check_angle.toDegrees())

anchors = AtomIndexSet()
anchors += residue.atom("H01")

residue.set(Angle("M03","O00","H01",ResNum(1)), 0.5*hoh, anchors)
#residue.set(Improper("M03","O00","H01","H02",ResNum(1)), SireMaths.Angle(0.0))

#check_angle = residue.measure("M03","O00","H01","H02")
#print "Improper angle = %f degrees." % check_angle.toDegrees()

print(tip4p)
print(residue)

PDB().write(tip4p, "tip4p.pdb")

