#!/bin/env python
# -*- coding: utf-8 -*-

from Sire.Mol import *   
from Sire.IO import *    
from Sire.Vol import *   
from Sire.FF import *    
from Sire.MM import *    
from Sire.Maths import * 
from Sire.Qt import *    
from Sire.Units import * 
from Sire.System import *
from Sire.Move import * 

t = QTime()

waters = PDB().read("test/io/water.pdb")

# only keep the first 50 water molecules
#waters2 = MoleculeGroup("waters")
#
#for i in range(0,50):
#    waters2.add( waters.moleculeAt(i) )
#
#waters = waters2

waters.setName("waters")

numbers = open("test/io/water.xsc", "r").readline().split()

for i in range(0,len(numbers)):
    numbers[i] = float(numbers[i])

mincoords = Vector( numbers[0], numbers[1], numbers[2] )
maxcoords = Vector( numbers[3], numbers[4], numbers[5] )

vol = PeriodicBox(mincoords, maxcoords)
vol = Cartesian()

t.start()
water = waters.moleculeAt(0).molecule() 

water = water.edit().atom( AtomName("O00") ) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom,  \
                                                   0.1550*kcal_per_mol)).molecule() \
                    .atom( AtomName("H01") ) \
                        .setProperty("charge", 0.520 * mod_electron).molecule() \
                    .atom( AtomName("H02") ) \
                        .setProperty("charge", 0.520 * mod_electron).molecule() \
                    .atom( AtomName("M03") ) \
                        .setProperty("charge", -1.04 * mod_electron).molecule() \
             .commit()

charges = water.property("charge")
ljs = water.property("LJ")

for i in range(0, waters.nMolecules()):
    water = waters.moleculeAt(i).molecule()

    water = water.edit().rename("T4P") \
                        .setProperty("charge", charges) \
                        .setProperty("LJ", ljs) \
                 .commit()

    waters.update(water)

ms = t.elapsed()
print("Parameterised all of the water molecules (in %d ms)!" % ms)

switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

cljff = InterCLJFF("CLJFF")

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

cljff.add( waters )

system = System()

system.add(cljff)
system.add(waters)

# constrain the identities of the first three waters
points = []

for i in range(0,3):
   points.append( waters.moleculeAt(i).molecule().evaluate().center() )

system.add( IdentityConstraint(points, waters) )

print("Calculating the starting energy... (should be -16364.5 kcal mol-1)")
print("...Initial energy = %s" % system.energy())

rbmc = RigidBodyMC(waters)
rbmc.setMaximumTranslation( 0.15*angstrom )
rbmc.setMaximumRotation( 15*degrees )
rbmc.setTemperature( 25*celsius )

#volmc = VolumeMove(waters)
#volmc.setMaximumVolumeChange( 0.1 * waters.nMolecules() * angstrom3 )
#volmc.setPressure( 1*atm )

moves = WeightedMoves()
#moves.add( volmc, 1 )
moves.add( rbmc, waters.nMolecules() )

# The simulation is now ready to run - create an output directory
# to collect the output
import os
import shutil

if os.path.exists("wateridentity_output"):
    print("Output directory already exists - removing it!")
    shutil.rmtree("wateridentity_output", True)

os.makedirs("wateridentity_output")

# Run 100 blocks of 10000 moves
for i in range(0,100):
    t.start()
    sim = Simulation.run(system, moves, 10000)
    ms = t.elapsed()

    print("Block %3d complete - took %d ms" % (i+1, ms))

    system = sim.system()
    moves = sim.moves()

    print("Block %3d: Energy = %f kcal mol-1" % (i+1, system.energy().to(kcal_per_mol)))

    PDB().write( system[MGName("waters")], "wateridentity_output/output%003d.pdb" % (i+1) )

    #xscfile = open("wateridentity_output/output%003d.xsc" % (i+1), "w")
    #vol = system.property("space")
    #maxcoords = vol.maxCoords()
    #
    #print >>xscfile,"%f %f %f  %f %f %f" % (mincoords[0], mincoords[1], mincoords[2], 
    #                                        maxcoords[0], maxcoords[1], maxcoords[2])

print("Simulation complete!")

