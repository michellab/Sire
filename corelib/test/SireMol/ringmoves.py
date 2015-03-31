
from Sire.Mol import *
from Sire.IO import *
from Sire.Units import *

import Sire.Stream

osel = Sire.Stream.load("test/io/osel.s3")
connectivity = osel.property("connectivity")

PDB().write(osel, "test0000.pdb")

map = { "weight function" : RelFromMass() }

bond = BondID( AtomName("C2"), AtomName("C9") )
angle = AngleID( AtomName("C2"), AtomName("C9"), AtomName("C8") )
dihedral = DihedralID( AtomName("C2"), AtomName("C9"), AtomName("C8"), AtomName("C5") )

def printPath(path):
    names = []

    for atom in path:
        names.append( osel.atom(atom).name().value() )

    return str(names)

def printPaths(paths):
    names = []

    for path in paths:
        names.append(printPath(path))

    return str(names)

# check that these are in the ring
print("%s is in a ring? %s { %s }" % (bond, connectivity.inRing(bond), \
                                        printPaths(connectivity.findPaths(bond.atom0(), bond.atom1()))))

print("%s is in a ring? %s { %s }" % (angle, connectivity.inRing(angle), \
                                        printPaths(connectivity.findPaths(angle.atom0(), angle.atom2()))))

print("%s is in a ring? %s { %s }" % (dihedral, connectivity.inRing(dihedral), \
                                        printPaths(connectivity.findPaths(dihedral.atom0(), dihedral.atom3()))))

osel = osel.move().change(bond, 0.3*angstrom, map).commit()

PDB().write(osel, "test0001.pdb")

osel = osel.move().change(bond, -0.3*angstrom, map).commit()

PDB().write(osel, "test0002.pdb")

osel = osel.move().change(angle, 15*degrees, map).commit()

PDB().write(osel, "test0003.pdb")

osel = osel.move().change(angle, -15*degrees, map).commit()

PDB().write(osel, "test0004.pdb")

osel = osel.move().change(dihedral, 30*degrees, map).commit()

PDB().write(osel, "test0005.pdb")

osel = osel.move().change(dihedral, -30*degrees, map).commit()

PDB().write(osel, "test0006.pdb")

