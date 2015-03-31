from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *
from Sire.IO import *
from Sire.Maths import *
from Sire.Units import *
from Sire.Qt import *

print("Loading a box of water...")
mols = PDB().read("test/io/water.pdb")

internalff = InternalFF()

print("Parameterising the molecules...")
for i in range(0, mols.nMolecules()):
    tip4p = mols.moleculeAt(i).molecule()

    bonds = Connectivity(tip4p)

    O00 = tip4p.atom( AtomName("O00") )
    H01 = tip4p.atom( AtomName("H01") )  
    H02 = tip4p.atom( AtomName("H02") )
    M03 = tip4p.atom( AtomName("M03") )

    tip4p = tip4p.edit().setProperty("connectivity", bonds).commit()

    bondfuncs = TwoAtomFunctions(tip4p)

    r = internalff.symbols().bond().r()

    bondfuncs.set( O00.index(), H01.index(), 300 * ( 0.9347 - r )**2 )
    bondfuncs.set( O00.index(), H02.index(), 300 * ( 0.9347 - r )**2 )

    anglefuncs = ThreeAtomFunctions(tip4p)

    theta = internalff.symbols().angle().theta()

    anglefuncs.set( H01.index(), O00.index(), H02.index(), 70 * ( (109.5*degrees).value() - theta )**2 )

    tip4p = tip4p.edit().setProperty( "bond", bondfuncs ) \
                        .setProperty( "angle", anglefuncs ) \
                 .commit()

    internalff.add( tip4p )

print("Calculating the intramolecular energy of the waters...")

t = QTime()
t.start()
print(internalff.energy())
ms = t.elapsed()

print(internalff.energy( internalff.components().bond() ))
print(internalff.energy( internalff.components().angle() ))

print("Calculation took %d ms" % ms)

print("\nManually calculating the energy for comparison...")

bndnrg = 0
angnrg = 0

t.start()
for i in range(0,mols.nMolecules()):
    tip4p = mols.moleculeAt(i).molecule()

    O00 = tip4p.atom(AtomName("O00")).property("coordinates")
    H01 = tip4p.atom(AtomName("H01")).property("coordinates")
    H02 = tip4p.atom(AtomName("H02")).property("coordinates")

    bndnrg += 300 * (0.9347 - Vector.distance(O00,H01))**2
    bndnrg += 300 * (0.9347 - Vector.distance(O00,H02))**2

    angnrg += 70 * (1.911135530933791 - Vector.angle(H01,O00,H02).value())**2

ms = t.elapsed()

print(bndnrg)
print(angnrg)
print(bndnrg + angnrg)

print("Manual calculation took %d ms" % ms)

print("\nTesting some moves...\n")

tip4p = internalff.molecule(MolIdx(0)).molecule()

new_tip4p = tip4p.move().change( BondID(AtomName("O00"),AtomName("H01")), 0.3 * angstrom ).commit()

internalff.update(new_tip4p)

t.start()
print(internalff.energy())
ms = t.elapsed()

print(internalff.energy( internalff.components().bond() ))
print(internalff.energy( internalff.components().angle() ))

print("Calculation took %d ms" % ms)

internalff.update(tip4p)

t.start()
print(internalff.energy())
ms = t.elapsed()

print(internalff.energy( internalff.components().bond() ))
print(internalff.energy( internalff.components().angle() ))

print("Calculation took %d ms" % ms)

print("\nTesting by running lots of moves...\n")

rand = RanGenerator()

t.start()

for i in range(0,1000):
    tip4p = internalff.molecule( MolIdx(rand.randInt(internalff.nMolecules()-1)) ).molecule()

    new_tip4p = tip4p.move().change( BondID(AtomName("O00"),AtomName("H01")), \
                                     rand.rand(-0.2, 0.2) * angstrom ) \
                            .change( AngleID(AtomName("H01"),AtomName("O00"),AtomName("H02")), \
                                     rand.rand(-5,5) * degrees ).commit()

    internalff.update(new_tip4p)

    #nrg = internalff.energy()

ms = t.elapsed()
print("\nMoves took %d ms\n" % ms)

t.start()
print(internalff.energy())
ms = t.elapsed()
print(internalff.energy( internalff.components().bond() ))
print(internalff.energy( internalff.components().angle() ))
print("Took %d ms" % ms)

internalff.mustNowRecalculateFromScratch()

print("\n", end=' ')
t.start()
print(internalff.energy())
ms = t.elapsed()
print(internalff.energy( internalff.components().bond() ))
print(internalff.energy( internalff.components().angle() ))
print("Took %d ms" % ms)
