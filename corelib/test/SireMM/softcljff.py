from Sire.MM import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.CAS import *
from Sire.Units import *

mol0 = Molecule("Sodium")
mol0 = mol0.edit().add(CGName("1")).add( AtomName("Na") ).molecule().commit()

mol1 = Molecule("Chloride")
mol1 = mol1.edit().add(CGName("1")).add( AtomName("Cl") ).molecule().commit()

mol0 = mol0.atom( AtomName("Na") ).edit() \
           .setProperty("coordinates", Vector(0,0,0)) \
           .setProperty("charge", 1 * mod_electron) \
           .setProperty("LJ", LJParameter(3.3284*angstrom, 0.0030*kcal_per_mol)) \
           .molecule().commit()

mol1 = mol1.atom( AtomName("Cl") ).edit() \
           .setProperty("coordinates", Vector(4,0,0)) \
           .setProperty("charge", -1 * mod_electron) \
           .setProperty("LJ", LJParameter(4.4010*angstrom, 0.1000*kcal_per_mol)) \
           .molecule().commit()

softclj = InterSoftCLJFF("SoftCLJ")

softclj.add( mol0 )
softclj.add( mol1 )

softclj.setAlpha(2, 0.0)
softclj.setAlpha(3, 0.5)
softclj.setAlpha(4, 1.0)
softclj.removeAlpha(0)

clj = InterCLJFF("CLJ")

clj.add( mol0 )
clj.add( mol1 )

def printComponents(energies):
    """Use this function to print out the energies that are in 'energies'"""
    
    components = energies.symbols()
    
    for component in components:
        print("%s == %f kcal mol-1" % (component, energies[component]))
    
    print("\n", end=' ')

printComponents( clj.energies() )

printComponents( softclj.energies() )
