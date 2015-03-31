
from Sire.Mol import *
from Sire.MM import *

c_1 = ( "C", 1 )
c_2 = ( "C", 2 )
h1_1 = ( "H1", 1 )
h2_1 = ( "H2", 1 )
h1_2 = ( "H1", 2 )
h2_2 = ( "H2", 2 )

bonds = MolBondInfo()

b = Bond( ("C",1), ("H",2) )
print(b)

bonds.addBond( (c_1,h1_1) )
bonds.addBond( (c_1,h2_1) )
bonds.addBond( (c_1,c_2) )
bonds.addBond( (c_2,h1_2) )
bonds.addBond( (c_2,h2_2) )

print(bonds.residuesBonded(1,2))

print(bonds.residuesBonded(1,1))

print(bonds.residuesBonded(1,3))

print(bonds.contains( (c_1,c_2) ))
print(bonds.contains( (h1_1,c_2) ))

it = bonds.begin()

print(it.isValid())
print(it.value())
print(it.key())

group_it = it.currentGroup()

print("Bonds in residues ",group_it.residueNumbers())
while group_it.isValid():
	print("%s = %s" % (group_it.value(), group_it.key()))
	next(group_it)

print("\nbonds")
print("%d bonds, %d intra-residue and %d inter-residue" % (bonds.nBonds(),bonds.nIntraBonds(),bonds.nInterBonds()))
print("intra-bonds")
it = bonds.intraBonds()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	it += 1

print("inter-bonds")
it = bonds.interBonds()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	next(it)

print("\nangles")

angs = MolAngleInfo()
angs.addAngle( (h1_1,c_1,c_2) )
angs.addAngle( (h2_1,c_1,c_2) )
angs.addAngle( (h1_1,c_1,h2_1) )
angs.addAngle( (c_1,c_2,h1_2) )
angs.addAngle( (h2_2,c_2,c_1) )
angs.addAngle( (h2_2,c_2,h1_2) )

print("%d angles, %d intra-residue and %d inter-residue" % (angs.nAngles(),angs.nIntraAngles(),angs.nInterAngles()))
print("(%d groups)" % angs.nGroups())

print("intra-angles")
it = angs.intraAngles()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	next(it)

print("inter-angles")
it = angs.interAngles()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	next(it)

print("\ndihedrals")

dihs = MolDihedralInfo()
dihs.addDihedral( (h2_1,c_1,c_2,h2_2) )
dihs.addDihedral( (h1_2,c_2,c_1,h1_1) )
dihs.addDihedral( (h1_2,c_2,c_1,h2_1) )
dihs.addDihedral( (h1_1,c_1,c_2,h2_2) )

print("%d dihedrals, %d intra-residue and %d inter-residue" % (dihs.nDihedrals(),dihs.nIntraDihedrals(),dihs.nInterDihedrals()))
print("(%d groups)" % dihs.nGroups())

print("intra-dihedrals")
it = dihs.intraDihedrals()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	next(it)

print("inter-dihedrals")
it = dihs.interDihedrals()
while (it.isValid()):
	print("%s = %s" % (it.value(),it.key()))
	next(it)

def getGroupComposition(it):
	
	if not it.isValid():
		print("WARNING: invalid iterator!")

	while(it.isValid()):
		print("Group %s : %s" % (it.value().groupID(),it.residueNumbers()))
		groupit = it.currentGroup()

		while groupit.isValid():
			print("%s = %s" % (groupit.value().index(), groupit.key()))
			next(groupit)

		it.nextGroup()

print("\nGroup composition for bonds...")
getGroupComposition( bonds.bonds() )

print("\nGroup composition for angles...")
getGroupComposition( angs.angles() )

print("\nGroup composition for dihedrals...")
getGroupComposition( dihs.dihedrals() )

print("Done!")




