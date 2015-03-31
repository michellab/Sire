
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *

import os

protomsdir = "%s/Work/ProtoMS" % os.getenv("HOME")

protoms = ProtoMS( "%s/protoms2" % protomsdir )

protoms.addParameterFile( "%s/parameter/amber99.ff" % protomsdir )
protoms.addParameterFile( "%s/parameter/gaff.ff" % protomsdir )
protoms.addParameterFile( "test/ff/sb2.ff" )

mol = PDB().readMolecule("test/io/sb2.pdb")

mol = mol.edit().rename("SB2").commit()

mol = protoms.parameterise(mol, ProtoMS.SOLUTE)

perturbations = mol.property("perturbations")

print(perturbations)

print(perturbations.requiredSymbols())
print(perturbations.requiredProperties())

lam = perturbations.symbols().Lambda()

system = System()

solute = MoleculeGroup("solute", mol)

system.add(solute)

system.setConstant(lam, 0.0)
system.add( PerturbationConstraint(solute) )

print(system.constraintsSatisfied())

for i in range(0,101,10):
    system.setConstant(lam, 0.01 * i)

    PDB().write(system.molecules(), "test_%003d.pdb" % i)

