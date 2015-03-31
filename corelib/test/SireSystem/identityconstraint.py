
from Sire.Base import *
from Sire.MM import *
from Sire.IO import *
from Sire.FF import *
from Sire.System import *
from Sire.Units import *
from Sire.Maths import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Qt import *

t = QTime()

cljff = InterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector( 18.3854,  18.66855,  18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("test/io/water.pdb")

print("Read in %d molecules!" % mols.nMolecules())

i = 0

t.start()
mol = mols.moleculeAt(0).molecule()

mol = mol.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.1550*kcal_per_mol)).molecule() \
                .atom( AtomName("H01") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("H02") ) \
                    .setProperty("charge", 0.520 * mod_electron).molecule() \
                .atom( AtomName("M03") ) \
                    .setProperty("charge", -1.04 * mod_electron).molecule() \
         .commit()

charges = mol.property("charge")
ljs = mol.property("LJ")

cljff.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = mol.edit().rename("T4P") \
                    .setProperty("charge", charges) \
                    .setProperty("LJ", ljs) \
             .commit()
   
    cljff.add(mol)    

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

def printMols( mols, molrange ):

   for i in molrange:
       mol_i = mols[ MolIdx(i) ]
       print(i, mol_i, mol_i.evaluate().center())

print("\nHERE ARE MOLECULES 0, 1, 10 and 20")
printMols( cljff, [0,1,10,20] )

system = System()

system.add(cljff)

t.start()
nrg = system.energy()
ms = t.elapsed()

print("\nEnergy = %f kcal mol-1 - took %d ms" % (nrg.to(kcal_per_mol), ms))

def testConstraint(points, system):

    print("Creating the identity constraint...")
    t.start()

    if len(points) == 0:
        idcons = IdentityConstraint( system[MGIdx(0)] )
    else:
        idcons = IdentityConstraint( points, system[MGIdx(0)] )

    ms = t.elapsed()
    print("...took %d ms" % ms)

    print("\nApplying the constraint...")
    t.start()

    mols = idcons.update(system)

    ms = t.elapsed()

    print("...took %d ms" % ms)

    return mols

centers = []

for i in range(0,system.nMolecules()):
    centers.append( cljff[MolIdx(i)].evaluate().center() )

def printMolecules(mols):
    for molnum in mols.molNums():
        mol = mols[molnum]
        print(molnum, mol, mol.evaluate().center())

print("\nHere are the coordinates of the centers of molecules 10 and 20")
print(centers[10], centers[20])

print("\nTesting centers[0], centers[1]")
mols = testConstraint( [centers[0], centers[1]], system )
print(mols)
printMolecules(mols)

print("\nTesting centers[10], centers[20]")
mols = testConstraint( [centers[10], centers[20]], system )
print(mols)
printMolecules(mols)

print("\nTesting centers[20], centers[10]")
mols = testConstraint( [centers[20], centers[10]], system )
print(mols)            
printMolecules(mols)

print("\nTesting []")
mols = testConstraint( [], system )
print(mols)
printMolecules(mols)

print("\nCOMPARISON TEST (centers[20], centers[10])")
idcons = IdentityConstraint( [centers[20], centers[10]], system[MGIdx(0)] )
mols = idcons.update(system)
print("\nDEFAULT")
print(mols)
printMolecules(mols)

print("\nFEWPOINTS")
idcons.useFewPointsAlgorithm()
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nMANYPOINTS")
idcons.useManyPointsAlgorithm()
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nREPEAT MANYPOINTS")
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nCOMPARISON TEST (centers[99])")
idcons = IdentityConstraint( [centers[99]], system[MGIdx(0)] )
mols = idcons.update(system)
print("\nDEFAULT")
print(mols)
printMolecules(mols)

print("\nSINGLEPOINT")
idcons.useSinglePointAlgorithm()
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nFEWPOINTS")
idcons.useFewPointsAlgorithm()
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nMANYPOINTS")
idcons.useManyPointsAlgorithm()
mols = idcons.update(system)
print(mols)
printMolecules(mols)

print("\nAPPLICATION TEST")

idcons = IdentityConstraint( [centers[99], centers[100], centers[101], centers[102]], system[MGIdx(0)] )
mols = idcons.update(system)
print("\nDEFAULT")
print(mols)
printMolecules(mols)

print("\nUPDATING SYSTEM")
print(system.version())
system.update(mols)
print(system.version())

mols = idcons.update(system)
print("\nPOST-UPDATE")
print(mols)
printMolecules(mols)
