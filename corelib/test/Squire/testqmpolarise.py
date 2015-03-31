
from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.Squire import *
from Sire.FF import *
from Sire.Vol import *
from Sire.Base import *
from Sire.Maths import *
from Sire.Units import *

import os

print("Loading the water molecules...")
waters = PDB().read("test/io/water.pdb")

tip4p = waters.moleculeAt(0).molecule()

protoms_dir = "%s/Work/ProtoMS" % os.getenv("HOME")

tip4p = tip4p.edit().rename("T4P").commit()

print("Parameterising the water as TIP4P...")
protoms = ProtoMS("%s/protoms2" % protoms_dir)
protoms.addParameterFile("%s/parameter/solvents.ff" % protoms_dir)

tip4p = protoms.parameterise(tip4p, ProtoMS.SOLVENT)

# set the polarisability of the atoms
tip4p = tip4p.edit() \
             .atom( AtomName("O00") ) \
                  .setProperty("polarisability", 0.465*angstrom3) \
                  .setProperty("fixed_charge", -0.656*mod_electron) \
                  .molecule() \
             .atom( AtomName("H01") ) \
                  .setProperty("polarisability", 0.135*angstrom3) \
                  .setProperty("fixed_charge", 0.328*mod_electron) \
                  .molecule() \
             .atom( AtomName("H02") ) \
                  .setProperty("polarisability", 0.135*angstrom3) \
                  .setProperty("fixed_charge", 0.328*mod_electron) \
                  .molecule() \
             .atom( AtomName("M03") ) \
                  .setProperty("polarisability", 0.465*angstrom3) \
                  .molecule() \
             .commit()

# Connect the atoms together...
connectivity = Connectivity(tip4p)
connectivity = connectivity.edit() \
                           .connect(AtomName("O00"), AtomName("H01")) \
                           .connect(AtomName("O00"), AtomName("H02")) \
                           .disconnect(AtomName("O00"), AtomName("M03")) \
                           .disconnect(AtomName("M03"), AtomName("H01")) \
                           .disconnect(AtomName("M03"), AtomName("H02")) \
                           .commit()

print(connectivity)

tip4p_chgs = tip4p.property("charge")
tip4p_fix = tip4p.property("fixed_charge")
tip4p_ljs = tip4p.property("LJ")
tip4p_pol = tip4p.property("polarisability")

for i in range(0, waters.nMolecules()):
    water = waters.moleculeAt(i).molecule()
    water = water.edit().setProperty("charge", tip4p_chgs) \
                        .setProperty("fixed_charge", tip4p_fix) \
                        .setProperty("LJ", tip4p_ljs) \
                        .setProperty("polarisability", tip4p_pol) \
                        .setProperty("connectivity", connectivity) \
                 .commit()

    waters.update(water)

print("Constructing the forcefields...")

switchfunc = HarmonicSwitchingFunction( 15*angstrom, 14.5*angstrom )
words = open("test/io/water.xsc", "r").readline().split()
space = PeriodicBox( Vector( float(words[0]), float(words[1]), float(words[2]) ),
                     Vector( float(words[3]), float(words[4]), float(words[5]) ) )

qm_water = waters.moleculeAt(0).molecule()
waters.remove(qm_water)

molpro = Molpro()

qmff = QMMMFF("qmmm")
qmff.setQuantumProgram(molpro)

qmff.add(qm_water, MGIdx(0))
qmff.add(waters, MGIdx(1))

system = System()
system.add(qmff)
system.setProperty("switchingFunction", switchfunc)
system.setProperty("space", space)

print("Calculating the system energies...")
print(system.energies())

polchgs = PolariseCharges(waters, qmff.components().total(),
                          CoulombProbe(1*mod_electron))

system.add(polchgs)
system.add(polchgs.selfEnergyFF())

print("Applying the polarisation constraint...")
system.applyConstraints()

print("Recalculating the energies...")
print(system.energies())

waters = system[qmff[MGIdx(1)].number()].molecules()

for molnum in waters.molNums():
    water = waters[molnum].molecule()
    print(molnum, water)
    print("MM charges\n",water.property("charge"), \
                         water.evaluate().charge({"charge":"charge"}))
    print("Fixed charges\n",water.property("fixed_charge"), \
                            water.evaluate().charge({"charge":"fixed_charge"}))
    print("Induced charges\n",water.property("induced_charge"), \
                              water.evaluate().charge({"charge":"induced_charge"}))
    print("New charges\n",water.property("charge"), \
                          water.evaluate().charge({"charge":"charge"}))

