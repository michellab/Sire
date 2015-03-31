
from Sire.IO import *
from Sire.Mol import *
from Sire.System import *
from Sire.MM import *
from Sire.FF import *
from Sire.Vol import *
from Sire.Base import *
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

pol_tip4p = waters.moleculeAt(0).molecule()
waters.remove(pol_tip4p)

cljff = InterGroupCLJFF("pol_tip4p-water")
cljff.add(pol_tip4p, MGIdx(0))
cljff.add(waters, MGIdx(1))

system = System()
system.add(cljff)

print(system.energies())

polchgs = PolariseCharges(cljff[MGIdx(0)], cljff.components().coulomb(),
                          CoulombProbe(1*mod_electron))

system.add(polchgs)
system.add(polchgs.selfEnergyFF())

print("Applying the polarisation constraint...")
system.applyConstraints()

print(system.energies())

pol_tip4p = system[MGIdx(0)][pol_tip4p.number()].molecule()

print("MM charges\n",tip4p.property("charge"), \
                     tip4p.evaluate().charge({"charge":"charge"}))
print("Fixed charges\n",pol_tip4p.property("fixed_charge"), \
                        pol_tip4p.evaluate().charge({"charge":"fixed_charge"}))
print("Induced charges\n",pol_tip4p.property("induced_charge"), \
                          pol_tip4p.evaluate().charge({"charge":"induced_charge"}))
print("New charges\n",pol_tip4p.property("charge"), \
                      pol_tip4p.evaluate().charge({"charge":"charge"}))

grid = RegularGrid(tip4p.evaluate().center(), 10, 1.0*angstrom)

def writePotential(molecule, filename):
    cljff = InterCLJFF()
    cljff.add(molecule)

    table = PotentialTable( cljff[MGIdx(0)], grid )

    cljff.potential(table, CoulombProbe(1*mod_electron))

    Cube(100*kcal_per_mol).write(table, cljff[MGIdx(0)], filename)

writePotential(tip4p, "tip4p.cube")
writePotential(pol_tip4p, "pol_tip4p.cube")

