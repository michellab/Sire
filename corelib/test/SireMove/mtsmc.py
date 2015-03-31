
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.FF import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *
from Sire.System import *
from Sire.Move import *
from Sire.Squire import *
from Sire.Base import *

import sys

timer = QTime()

#read in the solute molecule
print("Loading the tip4p solute...")
tip4p = PDB().read("test/io/tip4p.pdb")[0]

#read in all of the solvent molecules
print("Loading the solvent...")
solvent = PDB().read("test/io/water.pdb")

#renumber the the residues of the solvent so that they all have number 1
def renumberResidue(mols, newnum):
   
    newmols = []
    
    for mol in mols:
        editor = mol.edit()
        editor.setName("Water")
        editor.renumberResidue(mol.atoms()[0].resNum(), newnum)
        newmols.append( Molecule(editor) )

    return newmols

solvent = renumberResidue(solvent, ResNum(1))

#specify the space in which the molecules are placed
space = PeriodicBox( (-18.3854,-18.66855,-18.4445),
                     ( 18.3854, 18.66855, 18.4445) )

#specify the type of switching function to use
switchfunc = HarmonicSwitchingFunction(15.0, 14.5)

#create a forcefield to calculate the solvent-solvent energy
solventff = InterCLJFF(space, switchfunc)
solventff.setName("SolventFF")

#create a forcefield to calculate the MM solute-solvent energy
mmff = InterGroupCLJFF(space, switchfunc)
mmff.setName("MMFF")

#create a forcefield to calculate the QM solute-solvent energy
qmff = MolproFF(space, switchfunc)
qmff.setMolproExe( findExe("molpro") )
qmff.setName("QMFF")

qmff.setProgram("HF\nMP2")
qmff.setBasisSet("AVDZ")

#parametise each molecule and add it to the forcefields
print("Parametising the molecules...")

chgs = AtomicCharges( [0.0, 0.52 * mod_electron,
                            0.52 * mod_electron,
                           -1.04 * mod_electron] )

ljs = AtomicLJs( [ LJParameter( 3.15365 * angstrom, \
                                0.1550 * kcal_per_mol ), \
                   LJParameter.dummy(), \
                   LJParameter.dummy(), \
                   LJParameter.dummy() ] )

#first parametise the solute and add it to the right forcefields
tip4p.setProperty( "charges", chgs )
tip4p.setProperty( "ljs", ljs )

mmff.addToA( tip4p, {mmff.parameters().coulomb() : "charges",
                     mmff.parameters().lj() : "ljs"} )

qmff.addToQM( tip4p )

# QM energies are absolute, while MM energies are relative. We
# will set the QM energy relative to the QM energy of the gas-phase
# tip4p molecule
qmff.setEnergyOrigin( qmff.energy() )

#now parametise the solvent molecules and add them to the right
#forcefields
for mol in solvent:
    mol.setProperty( "charges", chgs )
    mol.setProperty( "ljs", ljs )

solventff.add( solvent, {solventff.parameters().coulomb() : "charges",
                         solventff.parameters().lj() : "ljs"} )
                         
mmff.addToB( solvent, {mmff.parameters().coulomb() : "charges",
                       mmff.parameters().lj() : "ljs"} )
                       
qmff.addToMM( solvent, {qmff.parameters().coulomb() : "charges"} )

#add all of the forcefields to a forcefields object
ffields = ForceFields()
ffields.add(solventff)
ffields.add(mmff)
ffields.add(qmff)

# there are two ways of calculating the total energy - either
# the slow way, using the QM solute-solvent interaction or
# the fast way, using the MM approximation of the solute-solvent
# interaction

e_fast = FFExpression( "e_fast", solventff.components().total() + 
                                 mmff.components().total() )

e_slow = FFExpression( "e_slow", solventff.components().total() + 
                                 mmff.components().lj() + 
                                 qmff.components().total() )

# we will perform a free energy simulation that perturbs from 
# e_slow to e_fast. This will allow us to calculate the correction
# free energy from e_fast to e_slow
#
#  We do this by using lambda to scale from e_slow to e_fast, e.g.
#   e_total = [(1-lambda) * e_slow] + [lambda * e_fast]

lam = Symbol("lambda")

e_total = FFExpression("e_total", ( (1-lam) * e_slow.function() ) 
                                    + ( lam * e_fast.function() ))

# the free energy gradient moving from e_slow to e_fast is given
# by the difference between e_fast and e_slow (dE / dlambda)
#
# Create an equation for this gradient so that we can monitor
# it during the simulation
de_by_dlam = FFExpression("de_by_dlam", e_fast.function() 
                                          - e_slow.function())

# add all of these energy components to the forcefields
ffields.add(e_fast)
ffields.add(e_slow)
ffields.add(e_total)
ffields.add(de_by_dlam)

# start at lambda = 0.0
ffields.setParameter(lam, 0.0)

# let's see what the initial values of these functions are...
#print "e_fast == %f kcal mol-1" % ffields.energy(e_fast)
#print "e_slow == %f kcal mol-1" % ffields.energy(e_slow)
#print "e_total == %f kcal mol-1" % ffields.energy(e_total)
#print "de_by_dlam == %f kcal mol-1" % ffields.energy(de_by_dlam)

#there are two groups of molecules, the solute and the solvent
# - lets create these groups (moves move random molecules from a group)
solvent_group = MoleculeGroup("solvent", solvent)
solute_group = MoleculeGroup("solute", tip4p)

all_mols = MoleculeGroup("all", solvent)
all_mols.add(tip4p)

groups = MoleculeGroups()
groups.add(solvent_group)
groups.add(solute_group)
groups.add(all_mols)

# create a monitor to monitor the average energies and RDFs
monitors = SystemMonitors()
monitors.set( e_slow.function(), MonitorEnergy(e_slow.function()) )
monitors.set( e_fast.function(), MonitorEnergy(e_fast.function()) )
monitors.set( de_by_dlam.function(), MonitorEnergy(de_by_dlam.function()) )
monitors.set( e_total.function(), MonitorEnergy(e_total.function()) )

rdfmonitor = RDFMonitor(tip4p, solvent_group)

oxygen = [("O00",ResNum(1))]
hydrogen = [("H01",ResNum(1)),("H02",ResNum(1))]

rdfmonitor.addRDF( oxygen, oxygen, RDF(0,15,150) )
rdfmonitor.addRDF( oxygen, hydrogen, RDF(0,15,150) )
rdfmonitor.addRDF( hydrogen, oxygen, RDF(0,15,150) )
rdfmonitor.addRDF( hydrogen, hydrogen, RDF(0,15,150) )

monitors.set( Symbol("RDF"), rdfmonitor )

# create a system that can be used to simulate these groups
# using the forcefields that we have constructed
system = System( groups, ffields, monitors )

# set the space in which the simulation will occur
system.setSpace(space)

# create a rigid body MC move that moves random solvent
# molecules, with the solvent molecules picked uniformly
mc = RigidBodyMC( PrefSampler(tip4p, 200.0, all_mols) )

# Create a multiple-time-step MC move that performs 500
# rigid body moves of the solvent using e_fast, then accepts
# or rejects the whole block of moves using e_total
mtsmc = MTSMC(mc, e_fast.function(), 500)

mtsmc.setEnergyComponent(e_total.function())

PDB().write(system.info().groups().molecules(), "test000.pdb")

for i in range(1,21):
    nmoves = 1000

    print("Block %d: Running %d moves" % (i, nmoves))

    timer.start()

    moves = system.run(mtsmc,nmoves)

    ms = timer.elapsed()

    print("%d moves took %d ms" % (nmoves, ms))

    PDB().write(system.info().groups().molecules(), "test%3.3d.pdb" % i)

mtsmc = moves.moves()[0].clone()

print("%d accepted, %d rejected, ratio == %f %%" % \
           (mtsmc.nAccepted(), mtsmc.nRejected(), mtsmc.acceptanceRatio()))

# check that the QM and MM energies have been conserved...
new_molpro = MolproFF(space, switchfunc)
new_molpro.setMolproExe( qmff.molproExe() )
new_molpro.setEnergyOrigin( qmff.energyOrigin() )

system_qmff = system.forceFields().forceField(qmff.ID())

new_molpro.addToQM( system_qmff.molecules(qmff.groups().qm()) )
new_molpro.addToMM( system_qmff.molecules(qmff.groups().mm()), 
                    {new_molpro.parameters().coulomb() : "charges"} )

print("QM energy = %f, new molproff = %f" % \
            ( system.forceFields().forceField(qmff.ID()).energy(),
              new_molpro.energy() ))

qmff.change(system.forceFields().molecules())

print("qmff == %f" % qmff.energy())

# get the RDFs
rdfmonitor = system.monitors().monitor(Symbol("RDF"))

rdf = rdfmonitor.getRDF( oxygen, oxygen )

print("OXYGEN-OXYGEN")

for point in rdf.normalise():
    print("%f  %f" % point)

print("\nOXYGEN-HYDROGEN")

rdf = rdfmonitor.getRDF( oxygen, hydrogen )

for point in rdf.normalise():
   print("%f  %f" % point)

print("\nHYDROGEN-OXYGEN")

rdf = rdfmonitor.getRDF( hydrogen, oxygen )

for point in rdf.normalise():
   print("%f  %f" % point)

print("\nHYDROGEN-HYDROGEN")

rdf = rdfmonitor.getRDF( hydrogen, hydrogen )

for point in rdf.normalise():
   print("%f  %f" % point)

