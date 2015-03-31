
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *

import copy

timer = QTime()

#read in all of the molecules
print("Loading the molecules...")
timer.start()
mols = PDB().read("test/io/water.pdb")

ms = timer.elapsed()
print("... took %d ms" % ms)

#specify the space in which the molecules are placed
space = Cartesian()

space = PeriodicBox( Vector(-18.3854,-18.66855,-18.4445),
                     Vector( 18.3854, 18.66855, 18.4445) )

#specify the type of switching function to use
switchfunc = HarmonicSwitchingFunction(80.0*angstrom)
switchfunc = HarmonicSwitchingFunction(15.0*angstrom, 14.5*angstrom)

#create a forcefield for the molecules
cljff = InterCLJFF(space, switchfunc)

coulff = InterCoulombFF(space, switchfunc)
ljff = InterLJFF(space, switchfunc)

cljff_a = InterCLJFF(space, switchfunc)
cljff_b = InterCLJFF(space, switchfunc)

cljff_a_b = InterGroupCLJFF(space, switchfunc)

ljff_a = InterLJFF(space, switchfunc)
ljff_b = InterLJFF(space, switchfunc)

ljff_a_b = InterGroupLJFF(space, switchfunc)

coulff_a = InterCoulombFF(space, switchfunc)
coulff_b = InterCoulombFF(space, switchfunc)

coulff_a_b = InterGroupCoulombFF(space, switchfunc)

#parametise each molecule and add it to the forcefield
print("Parametising the molecules...")

chgs = AtomicCharges( [0.0, 0.52 * mod_electron,
                            0.52 * mod_electron,
                           -1.04 * mod_electron] )

ljs = AtomicLJs( [ LJParameter( 3.15365 * angstrom, \
                                0.1550 * kcal_per_mol ), \
                   LJParameter.dummy(), \
                   LJParameter.dummy(), \
                   LJParameter.dummy() ] )

tip4p = False

timer.start()

rand = RanGenerator()

n_in_a = 0
n_in_b = 0

for mol in mols:
      
      mol.setProperty( "charges", chgs )
      mol.setProperty( "ljs", ljs )

      if (not tip4p):
          tip4p = mol

      #randomly divide the molecules into the two groups
      if (rand.rand() < 0.5):
          cljff_a.add(mol, {cljff_a.parameters().coulomb() : "charges",
                            cljff_a.parameters().lj() : "ljs"})
                            
          cljff_a_b.addTo(cljff_a_b.groups().A(),
                          mol, {cljff_a_b.parameters().coulomb() : "charges",
                                cljff_a_b.parameters().lj() : "ljs"})
      
          ljff_a.add(mol, {ljff_a.parameters().lj() : "ljs"})
          ljff_a_b.addTo(ljff_a_b.groups().A(),
                         mol, {ljff_a_b.parameters().lj() : "ljs"})

          coulff_a.add(mol, {coulff_a.parameters().coulomb() : "charges"})
          coulff_a_b.addTo(coulff_a_b.groups().A(),
                           mol, {coulff_a_b.parameters().coulomb() : "charges"})
                         
          n_in_a = n_in_a + 1
      else:
          cljff_b.add(mol, {cljff_b.parameters().coulomb() : "charges",
                            cljff_b.parameters().lj() : "ljs"})
                            
          cljff_a_b.addTo(cljff_a_b.groups().B(),
                          mol, {cljff_a_b.parameters().coulomb() : "charges",
                                cljff_a_b.parameters().lj() : "ljs"})
      
          ljff_b.add(mol, {ljff_b.parameters().lj() : "ljs"})
          ljff_a_b.addTo(ljff_a_b.groups().B(),
                         mol, {ljff_a_b.parameters().lj() : "ljs"})
      
          coulff_b.add(mol, {coulff_b.parameters().coulomb() : "charges"})
          coulff_a_b.addTo(coulff_a_b.groups().B(),
                           mol, {coulff_a_b.parameters().coulomb() : "charges"})

          n_in_b = n_in_b + 1

      cljff.add(mol, {cljff.parameters().coulomb() : "charges",
                      cljff.parameters().lj() : "ljs"})
                      
      coulff.add(mol, {coulff.parameters().coulomb() : "charges"})

      ljff.add(mol, {ljff.parameters().lj() : "ljs"})

ms = timer.elapsed()
print("... took %d ms" % ms)
      
print("(%d molecules in group A, %d in group B)" % (n_in_a, n_in_b))
      
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
nrg = coulff.energy()
ms = timer.elapsed()

print("InterCoulombFF ",coulff.energy(), "kcal mol-1")
print("   Coulomb = ", coulff.energy(coulff.components().coulomb()))

print("... took %d ms" % ms)


timer.start()
nrg = ljff.energy()
ms = timer.elapsed()

print("InterLJFF ",ljff.energy(), "kcal mol-1")
print("         LJ = ", ljff.energy(ljff.components().lj()))

print("... took %d ms" % ms)

print("Calculating the energy...")

# sum up the partial forcefields
timer.start()
nrg = cljff_a.energy() + cljff_b.energy() + cljff_a_b.energy()
ms = timer.elapsed()

print("CLJ_partials ",nrg,"kcal mol-1")
print("   Coulomb = ",cljff_a.energy(cljff_a.components().coulomb()) + \
                      cljff_b.energy(cljff_b.components().coulomb()) + \
                      cljff_a_b.energy(cljff_a_b.components().coulomb()))
                      
print("        LJ = ",cljff_a.energy(cljff_a.components().lj()) + \
                      cljff_b.energy(cljff_b.components().lj()) + \
                      cljff_a_b.energy(cljff_a_b.components().lj()))

print("... took %d ms" % ms)

timer.start()
nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()
ms = timer.elapsed()

print("LJ_partials ",nrg,"kcal mol-1")
print("... took %d ms" % ms)

timer.start()
nrg = coulff_a.energy() + coulff_b.energy() + coulff_a_b.energy()
ms = timer.elapsed()

print("Coulomb partials ",nrg,"kcal mol-1")
print("... took %d ms" % ms)

timer.start()

nmoves = 1000
for i in range(0,nmoves):
    cljff.change( tip4p )
    nrg = cljff.energy()

ms = timer.elapsed()

print("InterCLJFF ",cljff.energy(), "kcal mol-1")
print("... took %d ms (%f moves per second)" % (ms, nmoves*1000.0/ms))

tip4p = Molecule(tip4p.move().translate( (1.0,0.0,0.0) ))

timer.start()

cljff.change(tip4p)
ljff.change(tip4p)

ms = timer.elapsed()

print("Changing took %d ms" % ms)

timer.start()
nrg = cljff.energy()
ms = timer.elapsed()

print("InterCLJFF ",cljff.energy(), "kcal mol-1")
print("    Coulomb = ", cljff.energy(cljff.components().coulomb()))
print("         LJ = ", cljff.energy(cljff.components().lj()))

print("... took %d ms" % ms)

timer.start()
nrg = ljff.energy()
ms = timer.elapsed()

print("InterLJFF ",ljff.energy(), "kcal mol-1")
print("         LJ = ", ljff.energy(ljff.components().lj()))

print("... took %d ms" % ms)

timer.start()

nmoves = 1000
for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    cljff.change( tip4p )
    nrg = cljff.energy()

ms = timer.elapsed()

print("%d moves of InterCLJFF took %d ms" % (nmoves, ms))

timer.start()

nmoves = 1000

delta = 0

old_nrg = coulff.energy()
old_version = coulff.version()

for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    
    old_coulff = copy.copy(coulff)
    
    coulff.change( tip4p )
    nrg = coulff.energy()

    if (nrg == old_nrg):
       print("Energies are wrongly the same!!!")
    
    if (coulff.version() == old_version):
       print("Versions are wrongly the same!!!")
    
    coulff = old_coulff
    
    nrg = coulff.energy()
    
    if (nrg != old_nrg):
       print("Energies are wrongly different!!!")

    if (coulff.version() != old_version):
       print("Versions are wrongly different!!!")

ms = timer.elapsed()

print("%d moves of InterCoulombFF took %d ms" % (nmoves, ms))

timer.start()

nmoves = 1000

delta = 0

old_nrg = ljff.energy()
old_version = ljff.version()

for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    
    old_ljff = copy.copy(ljff)
    
    ljff.change( tip4p )
    nrg = ljff.energy()

    if (nrg == old_nrg):
       print("Energies are wrongly the same!!!")
    
    if (ljff.version() == old_version):
       print("Versions are wrongly the same!!!")
    
    ljff = old_ljff
    
    nrg = ljff.energy()
    
    if (nrg != old_nrg):
       print("Energies are wrongly different!!!")

    if (ljff.version() != old_version):
       print("Versions are wrongly different!!!")

ms = timer.elapsed()

print("%d moves of InterLJFF took %d ms" % (nmoves, ms))

timer.start()

nmoves = 1000

delta = 0

old_nrg = cljff_a.energy() + cljff_b.energy() + cljff_a_b.energy()

for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    
    old_cljff_a = copy.copy(cljff_a)
    old_cljff_b = copy.copy(cljff_b)
    old_cljff_a_b = copy.copy(cljff_a_b)
    
    cljff_a.change( tip4p )
    cljff_b.change( tip4p )
    cljff_a_b.change( tip4p )
    
    nrg = cljff_a.energy() + cljff_b.energy() + cljff_a_b.energy()

    if (nrg == old_nrg):
       print("Energies are wrongly the same!!!")
    
    cljff_a = old_cljff_a
    cljff_b = old_cljff_b
    cljff_a_b = old_cljff_a_b
    
    nrg = cljff_a.energy() + cljff_b.energy() + cljff_a_b.energy()
    
    if (nrg != old_nrg):
       print("Energies are wrongly different!!!")

ms = timer.elapsed()

print("%d moves of CLJ_partials took %d ms" % (nmoves, ms))

timer.start()

nmoves = 1000

delta = 0

old_nrg = coulff_a.energy() + coulff_b.energy() + coulff_a_b.energy()

for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    
    old_coulff_a = copy.copy(coulff_a)
    old_coulff_b = copy.copy(coulff_b)
    old_coulff_a_b = copy.copy(coulff_a_b)
    
    coulff_a.change( tip4p )
    coulff_b.change( tip4p )
    coulff_a_b.change( tip4p )
    
    nrg = coulff_a.energy() + coulff_b.energy() + coulff_a_b.energy()

    if (nrg == old_nrg):
       print("Energies are wrongly the same!!!")
    
    coulff_a = old_coulff_a
    coulff_b = old_coulff_b
    coulff_a_b = old_coulff_a_b
    
    nrg = coulff_a.energy() + coulff_b.energy() + coulff_a_b.energy()
    
    if (nrg != old_nrg):
       print("Energies are wrongly different!!!")

ms = timer.elapsed()

print("%d moves of Coulomb_partials took %d ms" % (nmoves, ms))

timer.start()

nmoves = 1000

delta = 0

old_nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()

for i in range(0,nmoves):
    tip4p = Molecule(tip4p.move().translate( (0.00001,0,0) ))
    
    old_ljff_a = copy.copy(ljff_a)
    old_ljff_b = copy.copy(ljff_b)
    old_ljff_a_b = copy.copy(ljff_a_b)
    
    ljff_a.change( tip4p )
    ljff_b.change( tip4p )
    ljff_a_b.change( tip4p )
    
    nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()

    if (nrg == old_nrg):
       print("Energies are wrongly the same!!!")
    
    ljff_a = old_ljff_a
    ljff_b = old_ljff_b
    ljff_a_b = old_ljff_a_b
    
    nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()
    
    if (nrg != old_nrg):
       print("Energies are wrongly different!!!")

ms = timer.elapsed()

print("%d moves of LJ_partials took %d ms" % (nmoves, ms))
