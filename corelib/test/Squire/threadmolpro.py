
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.FF import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Cluster import *
from Sire.Squire import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *

import time

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

molproexe = "../../../../../software/molpro/devel/molpro"

#create a forcefield for the molecules
molproff1 = MolproFF( space, switchfunc )
molproff2 = MolproFF( space, switchfunc )
molproff3 = MolproFF( space, switchfunc )

molproff1.setMolproExe(molproexe)
molproff2.setMolproExe(molproexe)
molproff3.setMolproExe(molproexe)

#parametise each molecule and add it to the forcefield
print("Parametising the molecules...")

chgs = AtomicCharges( [0.0, 0.52 * mod_electron, \
                            0.52 * mod_electron, \
                           -1.04 * mod_electron] )

ljs = AtomicLJs( [ LJParameter( 3.15365 * angstrom, \
                                0.1550 * kcal_per_mol ), \
                   LJParameter.dummy(), \
                   LJParameter.dummy(), \
                   LJParameter.dummy() ] )

timer.start()
for mol in mols:
      mol.setProperty( "charges", chgs )
      mol.setProperty( "ljs", ljs )

qm_mol = mols[0]
mm_mols = mols[1:]

molproff1.addToMM(mm_mols)
molproff2.addToMM(mm_mols)
molproff3.addToMM(mm_mols)

molproff1.addToQM(qm_mol)
molproff2.addToQM(qm_mol)
molproff3.addToQM(qm_mol)

ms = timer.elapsed()
print("... took %d ms" % ms)

timer.start()

#create a thread processor and calculate the energy in the background
threadproc1 = FFThreadProcessor(molproff1)
threadproc2 = FFThreadProcessor(molproff2)

active_threadproc1 = threadproc1.activate()
active_threadproc2 = threadproc2.activate()

print("Starting background calculation...")
active_threadproc1.recalculateEnergy()
active_threadproc2.recalculateEnergy()

print("Off it goes....")
print("Da de da da da...")

#create an FFProcessor, and place the cljff onto it...
ffproc1 = FFProcessor(molproff3)

print("Is active?", ffproc1.isActive())

active_ffproc1 = ffproc1.activate()

print("Is active?", ffproc1.isActive())

print("MAIN THREAD PROCESS")
print("Total energy == ",active_threadproc1.energy())
print("Total energy == ",active_threadproc2.energy())
print("Total energy == ",active_ffproc1.energy())

print("Took %d ms" % timer.elapsed())

