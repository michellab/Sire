
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.Base import *
from Sire.FF import *
from Sire.Maths import *
from Sire.CAS import *
from Sire.Units import *

import os

coords = "test/io/Bound_ligand1.pdb"
ff = "test/io/ligand1.template"
name = "dana"

coords = "test/io/ethane.pdb"
ff = "test/io/ethane.ff"
name = "ethane"

ethane = PDB().readMolecule(coords)

protoms_dir = "%s/Work/ProtoMS" % os.getenv("HOME")

protoms = ProtoMS("%s/protoms2" % protoms_dir)
protoms.addParameterFile("%s/parameter/amber99.ff" % protoms_dir)
protoms.addParameterFile("%s/parameter/gaff.ff" % protoms_dir)
protoms.addParameterFile(ff)

ethane = ethane.edit().rename(name).commit()
ethane = protoms.parameterise(ethane, ProtoMS.SOLUTE)

bonds = ethane.property("bond")
angles = ethane.property("angle")
dihedrals = ethane.property("dihedral")

print(bonds.potentials())
print(angles.potentials())
print(dihedrals.potentials())

#intraff = InternalFF("intraff")
intraclj = IntraCLJFF("intraclj")

#intraff.add(ethane)
intraclj.add(ethane)

solute = MoleculeGroup("solute", ethane)
solute.add(ethane)

system = System()
#system.add(intraff)
system.add(intraclj)
system.add(solute)

md = MolecularDynamics(solute, VelocityVerlet()) 
  #                     {"velocity generator" : MaxwellBoltzmann(25*celsius)})

md.setTimeStep(1*femtosecond)

PDB().write(system.molecules(), "test0000.pdb")

for i in range(1,250):
    md.move(system, 1)
    print(system.energy(), md.kineticEnergy(), (system.energy()+md.kineticEnergy()))
    PDB().write(system.molecules(), "test%0004d.pdb" % i)
