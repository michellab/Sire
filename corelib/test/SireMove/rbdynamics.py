
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
from Sire.Stream import *

import sys

t = QTime()

cljff = InterCLJFF("cljff")

c = Molecule("axes")
c = c.edit().add( CGName("1") ) \
            .add( AtomName("O") ).cutGroup() \
            .add( AtomName("+X") ).cutGroup() \
            .add( AtomName("-X") ).cutGroup() \
            .add( AtomName("+Y") ).cutGroup() \
            .add( AtomName("-Y") ).cutGroup() \
            .add( AtomName("+Z") ).cutGroup() \
            .add( AtomName("-Z") ).cutGroup() \
            .molecule().commit()

c = c.edit().atom( AtomName("O") ) \
            .setProperty("coordinates", Vector(0,0,0)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("+X") ) \
            .setProperty("coordinates", Vector(1,0,0)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("-X") ) \
            .setProperty("coordinates", Vector(-1,0,0)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("+Y") ) \
            .setProperty("coordinates", Vector(0,1,0)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("-Y") ) \
            .setProperty("coordinates", Vector(0,-1,0)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("+Z") ) \
            .setProperty("coordinates", Vector(0,0,1)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().atom( AtomName("-Z") ) \
            .setProperty("coordinates", Vector(0,0,-1)) \
            .setProperty("element", Element("Oxygen")) \
            .setProperty("charge", 1*mod_electron) \
            .setProperty("LJ", LJParameter(3*angstrom,0.5*kcal_per_mol)) \
            .molecule().commit()

c2 = c.edit().renumber().commit()

c2 = c2.edit().setProperty("charge", \
        AtomCharges( [ -1*mod_electron, -1*mod_electron, -1*mod_electron, \
          -1*mod_electron, -1*mod_electron, -1*mod_electron, \
          -1*mod_electron ] ) ).commit()

c2 = c2.move().translate( Vector(4,2,0) ).commit()

ion = Molecule("ion")
ion = ion.edit().add( CGName("1") ) \
                .add( AtomName("Na") ).molecule().commit()

salt = Molecule("salt")
salt = salt.edit().add( CGName("1") ) \
                  .add( AtomName("Na") ).cutGroup() \
                  .add( AtomName("Cl") ).molecule().commit()

ion = ion.edit().atom(AtomName("Na")) \
                .setProperty("coordinates", Vector(0,0,0)) \
                .setProperty("element", Element("Sodium")) \
                .setProperty("charge", 0.2*mod_electron) \
                .setProperty("LJ", LJParameter(3.0522*angstrom,0.4598*kcal_per_mol)) \
                .molecule().commit()

salt = salt.edit().atom(AtomName("Na")) \
                  .setProperty("coordinates", Vector(0,0,0)) \
                  .setProperty("element", Element("Sodium")) \
                  .setProperty("charge", 0.2*mod_electron) \
                  .setProperty("LJ", LJParameter(3.0522*angstrom,4.4598*kcal_per_mol)) \
                  .molecule() \
                  .atom(AtomName("Cl")) \
                  .setProperty("coordinates", Vector(0,2,0)) \
                  .setProperty("element", Element("Sodium")) \
                  .setProperty("charge", -0.2*mod_electron) \
                  .setProperty("LJ", LJParameter(3.0522*angstrom,4.4598*kcal_per_mol)) \
                  .molecule() \
                  .commit()

bonds = Connectivity(salt)
bonds = bonds.edit().connect( AtomName("Na"), AtomName("Cl") ).commit()

salt = salt.edit().setProperty("connectivity", bonds).commit()

bondfuncs = TwoAtomFunctions(salt)

internalff = InternalFF()

r = internalff.symbols().bond().r()         
bondfuncs.set( salt.atom(AtomName("Cl")).index(), salt.atom(AtomName("Na")).index(), 50 * ( 2.0 - r )**2 )

salt = salt.edit().setProperty("bond", bondfuncs).commit()

flexibility = Flexibility(salt)
flexibility.add( BondID(AtomIdx(0),AtomIdx(1)), 0.1*angstrom )
salt = salt.edit().setProperty("flexibility", flexibility).commit()

m0 = salt
m1 = salt.edit().renumber() \
         .atom(AtomName("Na")) \
         .setProperty("coordinates", Vector(4,2,1)) \
         .molecule().atom(AtomName("Cl")) \
         .setProperty("coordinates", Vector(4,4,1)) \
         .molecule().commit()

m2 = salt.edit().renumber() \
         .atom(AtomName("Na")) \
         .setProperty("coordinates", Vector(8,5,-1)) \
         .molecule().atom(AtomName("Cl")) \
         .setProperty("coordinates", Vector(8,7,-0.5)) \
         .molecule().commit()

m3 = salt.edit().renumber() \
         .atom(AtomName("Na")) \
         .setProperty("coordinates", Vector(12,0,0)) \
         .molecule().atom(AtomName("Cl")) \
         .setProperty("coordinates", Vector(14,0,0)) \
         .molecule().commit()

salt = MoleculeGroup("salt")
salt.add(m0)
salt.add(m1)
salt.add(m2)
salt.add(m3)
#salt.add(c)
#salt.add(c2)

cljff = InterCLJFF("salt-salt")
cljff.add(salt)
internalff.add(salt)

system = System()
system.add(salt)
system.add(cljff)
system.add(internalff)

system.setProperty("space", PeriodicBox( Vector(20,20,20) ) )

system.add( SpaceWrapper(Vector(0,0,0),salt) )

t.start()                                       
print("Initial energy = %s" % system.energy())
print("(took %d ms)" % t.elapsed())

mdmove = MolecularDynamics( salt, VelocityVerlet(), 
                            {"velocity generator":MaxwellBoltzmann(25*celsius)} )

mdmove = MolecularDynamics(salt, DLMRigidBody() )

mdmove.setTimeStep(1*femtosecond)
do_mc = False

intra_mcmove = InternalMove(salt)
inter_mcmove = RigidBodyMC(salt)

mcmove = WeightedMoves()
mcmove.add(intra_mcmove)
mcmove.add(inter_mcmove)

do_mc = False

hmcmove = HybridMC(salt, 4*femtosecond, 20)
do_hmc = True
do_hmc = False

print(system.property("space"))

print("\nMove 0")
print(system.energy())
print(mdmove.kineticEnergy())
print(system.energy() + mdmove.kineticEnergy())
PDB().write(system.molecules(), "test%0004d.pdb" % 0)

if do_mc:
    for i in range(1,1000):
        system = mcmove.move(system, 20, False)

        print(i, system.energy())

        PDB().write(system.molecules(), "test%0004d.pdb" % i)

elif do_hmc:
    for i in range(1,1000):
        hmcmove.move(system, 1)

        print(i, system.energy())
        PDB().write(system.molecules(), "test%0004d.pdb" % i)

else:
    for i in range(1,1000):
        print("\nmove %d" % (i))
        mdmove.move(system, 20)

        print(system.energy())
        print(mdmove.kineticEnergy())
        print(system.energy() + mdmove.kineticEnergy())

        PDB().write(system.molecules(), "test%0004d.pdb" % i)
