
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *

timer = QTime()

#read in all of the molecules
print("Loading the molecules...")
timer.start()
mols = PDB().read("test/io/water.pdb")

ms = timer.elapsed()
print("... took %d ms" % ms)

#specify the space in which the molecules are placed
space = Cartesian()

space = PeriodicBox(Vector(-18.3854,-18.66855,-18.4445), \
                    Vector( 18.3854, 18.66855, 18.4445))

#specify the type of switching function to use
switchfunc = HarmonicSwitchingFunction(80.0)
switchfunc = HarmonicSwitchingFunction(15.0, 14.5)

#create a forcefield for the molecules
cljff = InterCLJFF( space, switchfunc )

#parametise each molecule and add it to the forcefield
print("Parametising the molecules...")

chgs = AtomicCharges( [0.0, 0.52 * mod_electron,
                            0.52 * mod_electron,
                           -1.04 * mod_electron] )

ljs = AtomicLJs( [ LJParameter( 3.15365 * angstrom, 
                                0.1550 * kcal_per_mol ), 
                   LJParameter.dummy(), 
                   LJParameter.dummy(),
                   LJParameter.dummy() ] )

tip4p = False

timer.start()
for mol in mols:
      
      mol.setProperty( "charges", chgs )
      mol.setProperty( "ljs", ljs )

      if (not tip4p):
          tip4p = mol
      
      cljff.add(mol, {cljff.parameters().coulomb() : "charges",
                      cljff.parameters().lj() : "ljs"})

ms = timer.elapsed()
print("... took %d ms" % ms)
      
#now calculate the energy of the forcefield
print("Calculating the energy...")

timer.start()
nrg = cljff.energy()
ms = timer.elapsed()

print("InterCLJFF ",cljff.energy(), "kcal mol-1")
print("    Coulomb = ", cljff.energy(cljff.components().coulomb()))
print("         LJ = ", cljff.energy(cljff.components().lj()))

print("... took %d ms" % ms)

timer.start()

nmoves = 1000
for i in range(0,nmoves):
    cljff.change( tip4p )
    nrg = cljff.energy()

ms = timer.elapsed()

print("InterCLJFF ",cljff.energy(), "kcal mol-1")
print("... took %d ms (%f moves per second)" % (ms, nmoves*1000.0/ms))

