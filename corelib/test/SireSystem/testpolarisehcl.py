
from Sire.Mol import *
from Sire.MM import *
from Sire.FF import *
from Sire.System import *
from Sire.Maths import *
from Sire.Units import *

hcl = Molecule("HCl")
hcl = hcl.edit().add(CGName("1")).add(AtomName("H")) \
                .cutGroup().add(AtomName("Cl")).molecule().commit()

hcl = hcl.edit().setProperty("charge", AtomCharges([0.343*mod_electron,-0.343*mod_electron])) \
                .setProperty("element", AtomElements([Element("H"),Element("Cl")])) \
                .commit()

hcl = hcl.edit().atom(AtomName("H")).setProperty("polarisability", 10*angstrom3) \
         .molecule().atom(AtomName("Cl")).setProperty("polarisability", 10*angstrom3) \
         .molecule().commit()

hcl2 = hcl.edit().renumber().commit()

hcl = hcl.edit().setProperty("coordinates", AtomCoords([Vector(0),Vector(1.274,0,0)])).commit()
hcl2 = hcl2.edit().setProperty("coordinates", \
                       AtomCoords([Vector(2.5,0,0),Vector(3.774,0,0)])).commit()

hcl = hcl.edit().setProperty("connectivity", Connectivity(hcl)).commit()
hcl2 = hcl2.edit().setProperty("connectivity", Connectivity(hcl2)).commit()


hcl2_group = MoleculeGroup("HCl2", hcl2)

cff = InterCoulombFF("CoulombFF")
cff.add(hcl)
cff.add(hcl2)

print("Energy is %s : Should be about -10.9634 kcal mol-1" % cff.energy())

pottable = PotentialTable(hcl2_group)

cff.potential(pottable)

print("Potential is %s : Should be about [ -47.342 kcal mol-1, -15.380 kcal mol-1 ]" % \
          (str(pottable.getTable(hcl2.number()).toVector())))

pol = PolariseCharges(hcl2_group)

system = System()
system.add(cff)
system.add(hcl2_group)
system.add(pol.selfEnergyFF())
system.add(pol)

system.applyConstraints()

hcl2 = system[hcl2.number()].molecule()

print(system.energies())

print(hcl2.property("induced_charge"))
