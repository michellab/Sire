#!/bin/env python
# -*- coding: utf-8 -*-

###############################
#
# montecarlo.py
#
# This script performs a Monte Carlo
# simulation on a box of water.
#
# Rigid body moves are performed on
# the water molecules, which are chosen
# using a preferential sampling algorithm.
#
# An NPT ensemble is generated as 
# volume moves are attempted on the 
# periodic simulation box
#
#

# Import all of the necessary Sire modules
#################################################
from Sire.Mol import *      # Provides classes used to describe molecules
from Sire.IO import *       # Provides classes used to read and write molecules
from Sire.Vol import *      # Provides classes used to describe the periodic box
from Sire.FF import *       # Provides base classes used to describe the forcefields
from Sire.MM import *       # Provides classes used for molecular mechanics forcefields
from Sire.Maths import *    # Provides basic maths functions and classes
from Sire.Qt import *       # Exposes parts of the Qt library on which Sire is built
from Sire.Units import *    # Provides classes used to manipulate physical units (e.g. kcal mol-1)
from Sire.System import *   # Provides classes used to compose simulation systems
from Sire.Move import *     # Provides classes used to perform moves on the system

# Create a timer so we can measure how long everything takes
t = QTime()

# The input files are water.pdb, which contains the coordinates,
# and water.xsc, which contains the size of the periodic box

waters = PDB().read("water.pdb")   # PDB (from Sire.IO) reads PDB files

# The xsc file just contains 6 numbers, which are the minimum and maximum
# coordinates of the periodic box
numbers = open("water.xsc", "r").readline().split()

for i in range(0,len(numbers)):
    numbers[i] = float(numbers[i])

mincoords = Vector( numbers[0], numbers[1], numbers[2] )   # Vector (from Sire.Maths) describes
maxcoords = Vector( numbers[3], numbers[4], numbers[5] )   # a simple (x,y,z) vector

vol = PeriodicBox(mincoords, maxcoords)  # PeriodicBox (from Sire.Vol)
                                         # describes a periodic box

# Now lets parameterise the water molecules. We do this by editing one
# molecule, parameterising that, and then copying the parameters to
# the other water molecules
t.start()
water = waters.moleculeAt(0).molecule()   # .moleculeAt(0) returns the first molecule view.
                                          # .molecule() then converts this view to a plain Molecule
                                          # object (Molecule and MoleculeView are from Sire.Mol) 

# Molecules can be given arbitrary properties, with arbitrary names. To 
# assign a property, we need to .edit() the molecule. Once we have finished
# editing, we then .commit() the molecule.
#
# As we are editing atomic properties, we need to select each atom in turn
# (using .atom( AtomID ), where AtomID is something used to identify the atom,
#  in this case, AtomName, which is the name of the atom). We then set the
# property of the atom, via .setProperty( property_name, property_value ),
# before then returning to selecting the whole molecule  (via .molecule())
#
# The below lines of code perform the following sequence of events;
#  (1) editing of the water molecule is initiated, via .edit()
#  (2) the atom with name "O00" is selected, via .atom( AtomName("O00") )
#  (3) the "LJ" property of this atom is set to a LJ parameter,
#      via .setProperty("LJ", LJParameter(...))
#  (4) the whole molecule is then selected, via .molecule()
#  (5) the atom with name "H01" is selected, via .atom( AtomName("H01") )
#  (6) the "charge" property of this atom is set to 0.520 electron charges,
#      via .setProperty("charge", 0.520 * mod_electron)
#  (7) the whole molecule is then selected, via .molecule()
#  (8) the "H02" atom is selected, via .atom( AtomName("H02") )
#  (9) this is also then given a "charge" property of 0.520 electron charges
# (10) the whole molecule is then selected, via .molecule()
# (11) the "M03" atom is selected, via .atom( AtomName("M03") )
# (12) this atom is given a charge property of -1.04 electron charges
# (13) the whole molecule is selected, via .molecule()
# (14) finally, the editing is finished, via .commit(), which 
#      returns the, now parameterised, molecule 

# Note that this is really just one line of python - here it is broken
# up over multiple lines (using '\') for the purposes of readability

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

# The 'AtomName' class identifies atoms by name. It is provided by Sire.Mol, which
# also provides 'AtomNum' (atoms by number) and AtomIdx (atoms by index within the
# molecule, e.g. AtomIdx(2) would be the third atom in the molecule).
#
# The LJParameter class (provided by Sire.MM) provides a Lennard Jones parameter,
# with a sigma value (in angstrom) and an epsilon value (in kcal_per_mol)

# Now copy the "charge" and "LJ" properties from this water molecule
charges = water.property("charge")
ljs = water.property("LJ")

# Give these properties to all of the water molecules in 'waters' 
for i in range(0, waters.nMolecules()):
    water = waters.moleculeAt(i).molecule()

    # We set the properties, and also rename the molecule to "T4P" (via .rename("T4P"))
    water = water.edit().rename("T4P") \
                        .setProperty("charge", charges) \
                        .setProperty("LJ", ljs) \
                 .commit()

    # We have to update the 'waters' group with the new version of 
    # the molecule
    waters.update(water)

ms = t.elapsed()
print "Parameterised all of the water molecules (in %d ms)!" % ms

# We want to use a 15 A cutoff, which a 14.5 A feather - this is 
# provided by HarmonicSwitchingFunction (from Sire.MM), which takes
# the length of the cutoff and feather as arguments (angstrom is a length
# unit provided by Sire.Units)
switchfunc = HarmonicSwitchingFunction(15*angstrom, 14.5*angstrom)

# Now we need to create a forcefield that is used to calculate the
# intermolecular energy of the waters. The 'CLJ' forcefields 
# (Charge and Lennard Jones) can be used to calculate Coulomb
# and Lennard Jones energies. There are several CLJ forcefields,
# all provided in Sire.MM. The one we need here is called
# InterCLJFF, which is the ForceField to calculate the 
# Intermolecular Charge and Lennard Jones energy of all
# molecules contained within it.

cljff = InterCLJFF("CLJFF")

# We need to tell the forcefield which volume to use, and
# the size of the cutoff (the switching function)
cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

# We now need to add all of the water molecules to this forcefield,
# telling it to get the charges from the "charge" property, and the
# LJ parameters from the "LJ" property
cljff.add( waters, {"charge" : "charge", "LJ" : "LJ"} )

# Note that "charge" and "LJ" are the default names, so we could
# have just written 'cljff.add(waters)'

# Now create a simulation System - this holds all of the molecules,
# forcefields and monitors
system = System()

# Add the CLJ forcefield to the system
system.add(cljff)

# Add the molecule group containing the waters to the system - we will give
# this molecule group the name 'waters' so that we can extract it from
# the system at a later point in the script
waters.setName("waters")
system.add(waters)

# Add a constraint that waters are mapped back into the periodic box
system.add( SpaceWrapper(Vector(0,0,0), waters) )

print "Calculating the starting energy... (should be -16364.5 kcal mol-1)"
print "...Initial energy = %s" % system.energy()

# Now create the rigid body move (RigidBodyMC, from Sire.Move) that act on the waters
rbmc = RigidBodyMC(waters)

# Set the maximum amount by which to translate and rotate the water molecules
# to 0.25 A and 15 degrees (angstrom and degrees are from Sire.Units)
rbmc.setMaximumTranslation( 0.25*angstrom )
rbmc.setMaximumRotation( 15*degrees )

# Set the temperature of the moves to 25 Celsius (Celsius is from Sire.Units)
rbmc.setTemperature( 25*celsius )

# Now create the volume move (VolumeMove, from Sire.Move) that maintains a constant pressure
volmc = VolumeMove(waters)

# Set the maximum amount to change the volume to (nMolecules() / 10) A^3
volmc.setMaximumVolumeChange( 0.1 * waters.nMolecules() * angstrom3 )

# Set the pressure of the volume move to 1 atmosphere (atm is also from Sire.Units)
volmc.setPressure( 1*atm )

# Now create a WeightedMoves object (from Sire.Move) that will attempt
# one VolumeMove for every nMolecules() rigid body moves
moves = WeightedMoves()
moves.add( volmc, 1 )
moves.add( rbmc, 0.1 * waters.nMolecules() )

# The simulation is now ready to run - create an output directory
# to collect the output
import os
import shutil

if os.path.exists("montecarlo_output"):
    print "Output directory already exists - removing it!"
    shutil.rmtree("montecarlo_output", True)

os.makedirs("montecarlo_output")

# Run 200 blocks of 5000 moves
for i in range(0,200):
    t.start()
    sim = Simulation.run(system, moves, 5000)
    ms = t.elapsed()

    print "Block %3d complete - took %d ms" % (i+1, ms)

    # 'system' and 'moves' are the old copies - we need to
    # extract the new copies from 'sim'
    system = sim.system()
    moves = sim.moves()

    # Print out the energy of the system - note that you must explicitly
    # state the units you want the energy to be output it, using .to(), 
    # in this case we are using kcal mol-1, so we have .to(kcal_per_mol)
    print "Block %3d: Energy = %f kcal mol-1" % (i+1, system.energy().to(kcal_per_mol))

    # Now print out a new PDB of the waters - note that MGName("waters") selects
    # the molecule group with name 'waters' 
    PDB().write( system[MGName("waters")], "montecarlo_output/output%003d.pdb" % (i+1) )

    # Now write out the new XSC file for the simulation box
    xscfile = open("montecarlo_output/output%003d.xsc" % (i+1), "w")
 
    # the volume is held in the "space" property of the system
    vol = system.property("space")

    print "Simulation box: %s" % vol

    mincoords = vol.minCoords()
    maxcoords = vol.maxCoords()

    print >>xscfile,"%f %f %f  %f %f %f" % (mincoords[0], mincoords[1], mincoords[2], 
                                            maxcoords[0], maxcoords[1], maxcoords[2])

    xscfile.close()

print "Simulation complete!"

